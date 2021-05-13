#if __APPLE__
#include <GLUT/glut.h>
#else
#include <windows.h>
#include <GL/glut.h>
#endif

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <algorithm>
using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// References:
// https://github.com/yuanming-hu/taichi_mpm/blob/master/mls-mpm88-explained.cpp
// https://lucasschuermann.com/writing/implementing-sph-in-2d for visualization

// Particle representation
struct Particle
{
	Particle(double _x, double _y, float r, float g, float b) : x(_x, _y), v(0.0, 0.0), F(Matrix2d::Identity()), C(Matrix2d::Zero()), Jp(1.0), color(r, g, b) {}
	Vector2d x, v; //position and velocity
	Matrix2d F, C; //deformation gradient, APIC momentum
	double Jp; //determinant of deformation gradient, which is volume
	Vector3f color;
};

// Grid representation
/*
struct Cell
{
	Cell() : v(0.f, 0.f), mass(0.f) {}
	Vector2f v;
	double mass;

};
*/

struct Rect
{
	Rect(float _x1, float _y1, float _x2, float _y2) : x1(_x1), y1(_y1), x2(_x2), y2(_y2) {}
	float x1, y1, x2, y2;
};

struct Tri
{
	Tri(double ax, double ay, double bx, double by, double cx, double cy) : a(ax, ay), b(bx, by), c(cx, cy) {}
	Vector2d a, b, c;
};

// Granularity
const static int MAX_PARTICLES = 25000;
const static int BLOCK_PARTICLES = 1000;		// number of particles added in a block
int NUM_PARTICLES = 0;					// keeps track of current number of particles
const static int GRID_RES = 80;				// grid dim of one side
const static int NUM_CELLS = GRID_RES * GRID_RES;	// number of cells in the grid
const static double DT = 0.0001;			// integration timestep
const static double DX = 1.0 / GRID_RES;
const static double INV_DX = 1.0 / DX;
const static double EPS = 0.005;

// Data structures
static vector<Particle> particles;
// Vector3: [velocity_x, velocity_y, mass]
static Vector3d grid[GRID_RES + 1][GRID_RES + 1];
//static Cell grid[GRID_RES][GRID_RES];
static vector<Rect> rectangles;
static vector<Tri> triangles;

// Simulation params
const static double MASS = 1.0;					// mass of one particle
const static double VOL = 1.0;					// volume of one particle
const static double HARD = 10.0;					// snow hardening factor, typically between 3-10
const static double E = 200;				// Young's Modulus, resistance to fracture
const static double NU = 0.0;					// Poisson ratio
const static double MU_V = 0.5;					// friction damping

// Initial Lame params
const static double MU_0 = E / (2 * (1 + NU));
const static double LAMBDA_0 = (E * NU) / ((1 + NU) * (1 - 2 * NU));

// Render params
const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
//const static double VIEW_WIDTH = 1.5 * 800;
//const static double VIEW_HEIGHT = 1.5 * 600;

// Image output params
static int step = 0;	// current simulation step
const static int FPS = 300; // fps of output video
const static int INV_FPS = (1.0 / DT) / FPS;

static int maxPartsPerGrid = 0;

static bool tetris = false; //turn on tetris mode, just for fun

// yay add more particles randomly in a square
void addParticles(double xcenter, double ycenter)
{
	for (int i = 0; i < BLOCK_PARTICLES; i++) {
		particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.08 + xcenter, (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.08 + ycenter, 0.3f, 0.8f, 0.3f));
	}
	NUM_PARTICLES += BLOCK_PARTICLES;
}

// a friend coerced me into adding this "feature"
void addTetris(double xcenter, double ycenter)
{
	int type = rand() % 6;
	for (int i = 0; i < BLOCK_PARTICLES / 4; i++) {
		switch (type) {
		case 0:
			//i block
			for (int j = 0; j < 4; j++) {
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * j), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter, 0.0f, 1.0f, 1.0f));
			}
			break;
		case 1:
			//l block
			for (int j = 0; j < 3; j++) {
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * j), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter, 1.0f, 0.5f, 0.0f));
			}
			particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * 2), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter + (0.06 * 1), 1.0f, 0.5f, 0.0f));
			break;
		case 2:
			//j block
			for (int j = 0; j < 3; j++) {
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * j), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter, 0.0f, 0.0f, 1.0f));
			}
			particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter, (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter + (0.06 * 1), 0.0f, 0.0f, 1.0f));
			break;
		case 3:
			//t block
			for (int j = 0; j < 2; j++) {
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * j), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter, 0.5f, 0.0f, 1.0f));
			}
			particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * 2), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter + (0.06 * 1), 0.5f, 0.0f, 1.0f));
			break;
		case 4:
			//o block
			for (int j = 0; j < 2; j++) {
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * j), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter, 1.0f, 1.0f, 0.0f));
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * j), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter + (0.06 * 1), 1.0f, 1.0f, 0.0f));
			}
			break;
		case 5:
			//s block
			for (int j = 0; j < 2; j++) {
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * j), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter, 0.0f, 1.0f, 0.0f));
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * (j + 1)), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter + (0.06 * 1), 0.0f, 1.0f, 0.0f));
			}
			break;
		case 6:
			//z block
			for (int j = 0; j < 2; j++) {
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * (j + 1)), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter, 1.0f, 0.0f, 0.0f));
				particles.push_back(Particle((((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + xcenter + (0.06 * j), (((double)rand() / (double)RAND_MAX) * 2 - 1) * 0.03 + ycenter + (0.06 * 1), 1.0f, 0.0f, 0.0f));
			}
			break;
		default:
			cout << "Error in adding tetris\n";
			break;
		}
	}
}

void InitMPM(void)
{
	cout << "initializing mpm with " << BLOCK_PARTICLES << " particles" << endl;
	addParticles(0.55, 0.45);
	addParticles(0.45, 0.65);
	addParticles(0.55, 0.85);
}

void addSlope()
{
	// Add a ramp in the lower left corner
	triangles.push_back(Tri(0.05, 0.05, 0.05, 0.5, 0.6, 0.05));
}

//Initialize obstacles. Makes a fun shape.
void addRectangles()
{
	rectangles.push_back(Rect(0.1f, 0.2f, 0.14f, 0.5f));

	rectangles.push_back(Rect(0.18f, 0.46f, 0.30f, 0.5f));
	rectangles.push_back(Rect(0.18f, 0.37f, 0.22f, 0.46f));
	rectangles.push_back(Rect(0.18f, 0.33f, 0.3f, 0.37f));
	rectangles.push_back(Rect(0.26f, 0.24f, 0.3f, 0.33f));
	rectangles.push_back(Rect(0.18f, 0.2f, 0.3f, 0.24f));

	rectangles.push_back(Rect(0.34f, 0.33f, 0.38f, 0.5f));
	rectangles.push_back(Rect(0.38f, 0.33f, 0.42f, 0.37f));
	rectangles.push_back(Rect(0.42f, 0.2f, 0.46f, 0.5f));

	rectangles.push_back(Rect(0.5f, 0.46f, 0.62f, 0.5f));
	rectangles.push_back(Rect(0.5f, 0.2f, 0.54f, 0.46f));
	rectangles.push_back(Rect(0.54f, 0.33f, 0.62f, 0.37f));
	rectangles.push_back(Rect(0.58f, 0.24f, 0.62f, 0.33f));
	rectangles.push_back(Rect(0.54f, 0.2f, 0.62f, 0.24f));

	rectangles.push_back(Rect(0.66f, 0.33f, 0.7f, 0.5f));
	rectangles.push_back(Rect(0.7f, 0.33f, 0.74f, 0.37f));
	rectangles.push_back(Rect(0.74f, 0.2f, 0.78f, 0.5f));
}

int SideTest(Vector2d x, Vector2d A, Vector2d B) {
	double tmp = (B.x() - A.x()) * (x.y() - A.y()) - (B.y() - A.y()) * (x.x() - A.x());
	return (tmp > 0) ? 1 : ((tmp < 0) ? -1 : 0);
}

// checks for collisions with any object and modifies velocity
void CheckCollisions(double x, double y, Vector3d& g) {
	Vector2d res(0.0, 0.0);
	for (Rect r : rectangles) {
		if (y > r.y1 && y < r.y2 && x < r.x2 + EPS && x > r.x1 - EPS) {
			g[0] = 0.0;
			return;
		}
		if (x > r.x1 && x < r.x2 && y < r.y2 + EPS && y > r.y1) {
			g[1] = 0.0;
			return;
		}
	}
	//if I had multiple triangles I would be testing all sides but since I know it can only collide with BC here I only coded BC
	for (Tri t : triangles) {
		Vector2d pos(x, y);
		Vector2d vec1 = pos - t.b;
		Vector2d bc = t.c - t.b;
		double len = bc.dot(vec1) / bc.squaredNorm();
		double dist = (vec1 - (len * bc)).norm();
		if (len >= 0 && len <= 1 && dist < EPS*2) {
			//Vector2d vtmp(g[0], g[1]);
			//res = vtmp - len * (t.c - t.b);
			res = len * bc;
			// slippery
			g[0] = res.x();
			g[1] = res.y();
			return;
		}
	}
	return;
}

void P2G(void)
{
	memset(grid, 0, sizeof(grid));
	for (Particle& p : particles) {
		Vector2i base_coord = (p.x * INV_DX - Vector2d(0.5, 0.5)).cast<int>();
		Vector2d fx = p.x * INV_DX - base_coord.cast<double>();

		// Quadratic B-spline weight kernels [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Eqn. 123, with x=fx, fx-1,fx-2]
		// This is referred to as N in the MLS-MPM paper
		Vector2d onehalf(1.5, 1.5); //have to make these so eigen doesn't give me errors
		Vector2d one(1.0, 1.0);
		Vector2d half(0.5, 0.5);
		Vector2d threequarters(0.75, 0.75);
		Vector2d tmpa = onehalf - fx;
		Vector2d tmpb = fx - one;
		Vector2d tmpc = fx - half;
		Vector2d w[3] = {
		  0.5 * (tmpa.cwiseProduct(tmpa)),
		  threequarters - (tmpb.cwiseProduct(tmpb)),
		  0.5 * (tmpc.cwiseProduct(tmpc))
		};

		// Snow
		// Compute current Lamé parameters [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Eqn. 87]
		double e = exp(HARD * (1.0 - p.Jp));
		double mu = MU_0 * e;
		double lambda = LAMBDA_0 * e;

		// Current volume
		double J = p.F.determinant();

		// Polar decomposition for fixed corotated model, https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf paragraph after Eqn. 45
		JacobiSVD<Matrix2d> svd(p.F, ComputeFullU | ComputeFullV);
		Matrix2d r = svd.matrixU() * svd.matrixV().transpose();
		//Matrix2d s = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();

		// [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Paragraph after Eqn. 176]
		double Dinv = 4 * INV_DX * INV_DX;
		// [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Eqn. 52]
		Matrix2d PF_0 = (2 * mu) * (p.F - r) * p.F.transpose();
		double pf1tmp = (lambda * (J - 1) * J);
		Matrix2d PF_1;
		PF_1 << pf1tmp, pf1tmp,
			pf1tmp, pf1tmp;
		Matrix2d PF = PF_0 + PF_1;

		// Cauchy stress times dt and inv_dx
		Matrix2d stress = -(DT * VOL) * (Dinv * PF);

		// Fused APIC momentum + MLS-MPM stress contribution
		// See https://yzhu.io/publication/mpmmls2018siggraph/paper.pdf
		// Eqn 29
		Matrix2d affine = stress + MASS * p.C;

		/*
		//Neo Hookean
		// [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Paragraph after Eqn. 176]
		double Dinv = 4 * INV_DX * INV_DX;
		// https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf equation 48
		Matrix2d P_0 = mu * (p.F - p.F.transpose().inverse());
		Matrix2d P_1 = lambda * log(J) * p.F.transpose().inverse();
		Matrix2d P = P_0 + P_1;

		// equation 38, MPM course
		Matrix2d stress = (1.0f / J) * (P * p.F.transpose());

		// Cauchy stress times dt and inv_dx
		stress = -(DT * VOL) * (Dinv * stress);

		// Fused APIC momentum + MLS-MPM stress contribution
		// See https://yzhu.io/publication/mpmmls2018siggraph/paper.pdf
		// Eqn 29
		Matrix2d affine = stress + MASS * p.C;
		*/

		// P2G
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Vector2d dpos = (Vector2d(i, j) - fx) * DX;
				// Translational momentum
				Vector3d mass_x_velocity(p.v.x() * MASS, p.v.y() * MASS, MASS);
				Vector2d tmp = affine * dpos;
				grid[base_coord.x() + i][base_coord.y() + j] += (
					w[i].x() * w[j].y() * (mass_x_velocity + Vector3d(tmp.x(), tmp.y(), 0))
					);
			}
		}
	}
}

void P2GJello(void)
{
	memset(grid, 0, sizeof(grid));
	for (Particle& p : particles) {
		Vector2i base_coord = (p.x * INV_DX - Vector2d(0.5, 0.5)).cast<int>();
		Vector2d fx = p.x * INV_DX - base_coord.cast<double>();

		// Quadratic B-spline weight kernels [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Eqn. 123, with x=fx, fx-1,fx-2]
		// This is referred to as N in the MLS-MPM paper
		Vector2d onehalf(1.5, 1.5); //have to make these so eigen doesn't give me errors
		Vector2d one(1.0, 1.0);
		Vector2d half(0.5, 0.5);
		Vector2d threequarters(0.75, 0.75);
		Vector2d tmpa = onehalf - fx;
		Vector2d tmpb = fx - one;
		Vector2d tmpc = fx - half;
		Vector2d w[3] = {
		  0.5 * (tmpa.cwiseProduct(tmpa)),
		  threequarters - (tmpb.cwiseProduct(tmpb)),
		  0.5 * (tmpc.cwiseProduct(tmpc))
		};

		// Current volume
		double J = p.F.determinant();
		// [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Paragraph after Eqn. 176]
		double Dinv = 4 * INV_DX * INV_DX;
		// https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf equation 48
		Matrix2d P_0 = MU_0 * (p.F - p.F.transpose().inverse());
		Matrix2d P_1 = LAMBDA_0 * log(J) * p.F.transpose().inverse();
		Matrix2d P = P_0 + P_1;

		// equation 38, MPM course
		Matrix2d stress = (1.0f / J) * (P * p.F.transpose());

		// Cauchy stress times dt and inv_dx
		stress = -(DT * VOL) * (Dinv * stress);

		// Fused APIC momentum + MLS-MPM stress contribution
		// See https://yzhu.io/publication/mpmmls2018siggraph/paper.pdf
		// Eqn 29
		Matrix2d affine = stress + MASS * p.C;

		// P2G
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Vector2d dpos = (Vector2d(i, j) - fx) * DX;
				// Translational momentum
				Vector3d mass_x_velocity(p.v.x() * MASS, p.v.y() * MASS, MASS);
				Vector2d tmp = affine * dpos;
				grid[base_coord.x() + i][base_coord.y() + j] += (
					w[i].x() * w[j].y() * (mass_x_velocity + Vector3d(tmp.x(), tmp.y(), 0))
					);
			}
		}
	}
}

void UpdateGridVelocity(void) {
	// For all grid nodes
	for (int i = 0; i <= GRID_RES; i++) {
		for (int j = 0; j <= GRID_RES; j++) {
			auto& g = grid[i][j];
			// No need for epsilon here
			if (g[2] > 0) {
				if (g[2] > maxPartsPerGrid) {
					maxPartsPerGrid = g[2];
				}
				// Normalize by mass
				g /= g[2];
				// Gravity
				g += DT * Vector3d(0, -200, 0);

				// boundary thickness
				double boundary = 0.05;
				// Node coordinates
				double x = (double)i / GRID_RES;
				double y = (double)j / GRID_RES;

				// Sticky boundary
				if (x < boundary || x > 1 - boundary) {
					g[0] = 0.0;
				}
				// Separate boundary
				if (y < boundary || y > 1 - boundary) {
					g[0] *= MU_V;
					g[1] = 0.0;
				}

				CheckCollisions(x, y, g);
			}
		}
	}
}

void G2P(void)
{
	for (Particle& p : particles) {
		// element-wise floor
		Vector2i base_coord = (p.x * INV_DX - Vector2d(0.5, 0.5)).cast<int>();
		Vector2d fx = p.x * INV_DX - base_coord.cast<double>();

		// Quadratic kernels [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Eqn. 123, with x=fx, fx-1,fx-2]
		Vector2d onehalf(1.5, 1.5); //have to make these so eigen doesn't give me errors
		Vector2d one(1.0, 1.0);
		Vector2d half(0.5, 0.5);
		Vector2d threequarters(0.75, 0.75);
		Vector2d tmpa = onehalf - fx;
		Vector2d tmpb = fx - one;
		Vector2d tmpc = fx - half;
		Vector2d w[3] = {
		  0.5 * (tmpa.cwiseProduct(tmpa)),
		  threequarters - (tmpb.cwiseProduct(tmpb)),
		  0.5 * (tmpc.cwiseProduct(tmpc))
		};

		p.C = Matrix2d::Zero();
		p.v = Vector2d::Zero();

		// constructing affine per-particle momentum matrix from APIC / MLS-MPM.
		// see APIC paper (https://www.math.ucla.edu/~jteran/papers/JSSTS15.pdf), page 6
		// below equation 11 for clarification. this is calculating C = B * (D^-1) for APIC equation 8,
		// where B is calculated in the inner loop at (D^-1) = 4 is a constant when using quadratic interpolation functions
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Vector2d dpos = (Vector2d(i, j) - fx);
				Vector3d curr = grid[base_coord.x() + i][base_coord.y() + j];
				Vector2d grid_v(curr.x(), curr.y());
				double weight = w[i].x() * w[j].y();
				// Velocity
				p.v += weight * grid_v;
				// APIC C, outer product of weighted velocity and dist, paper equation 10
				p.C += 4 * INV_DX * ((weight * grid_v) * dpos.transpose());
			}
		}

		// Advection
		//double tempy = p.x.y();
		p.x += DT * p.v;
		//cout << "change" << tempy - p.x.y() << endl;
		//cout << "velocity " <<  p.v << endl;

		// MLS-MPM F-update eqn 17
		Matrix2d F = (Matrix2d::Identity() + DT * p.C) * p.F;
		JacobiSVD<Matrix2d> svd(F, ComputeFullU | ComputeFullV);
		Matrix2d svd_u = svd.matrixU();
		Matrix2d svd_v = svd.matrixV();
		// Snow Plasticity
		Vector2d sigvalues = svd.singularValues().array().min(1.0f + 7.5e-3).max(1.0 - 2.5e-2);
		Matrix2d sig = sigvalues.asDiagonal();

		double oldJ = F.determinant();
		F = svd_u * sig * svd_v.transpose();

		double Jp_new = min(max(p.Jp * oldJ / F.determinant(), 0.6), 20.0);
		p.Jp = Jp_new;
		p.F = F;
	}
}

void G2PJello(void)
{
	for (Particle& p : particles) {
		// element-wise floor
		Vector2i base_coord = (p.x * INV_DX - Vector2d(0.5, 0.5)).cast<int>();
		Vector2d fx = p.x * INV_DX - base_coord.cast<double>();

		// Quadratic kernels [https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Eqn. 123, with x=fx, fx-1,fx-2]
		Vector2d onehalf(1.5, 1.5); //have to make these so eigen doesn't give me errors
		Vector2d one(1.0, 1.0);
		Vector2d half(0.5, 0.5);
		Vector2d threequarters(0.75, 0.75);
		Vector2d tmpa = onehalf - fx;
		Vector2d tmpb = fx - one;
		Vector2d tmpc = fx - half;
		Vector2d w[3] = {
		  0.5 * (tmpa.cwiseProduct(tmpa)),
		  threequarters - (tmpb.cwiseProduct(tmpb)),
		  0.5 * (tmpc.cwiseProduct(tmpc))
		};

		p.C = Matrix2d::Zero();
		p.v = Vector2d::Zero();

		// constructing affine per-particle momentum matrix from APIC / MLS-MPM.
		// see APIC paper (https://www.math.ucla.edu/~jteran/papers/JSSTS15.pdf), page 6
		// below equation 11 for clarification. this is calculating C = B * (D^-1) for APIC equation 8,
		// where B is calculated in the inner loop at (D^-1) = 4 is a constant when using quadratic interpolation functions
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Vector2d dpos = (Vector2d(i, j) - fx);
				Vector3d curr = grid[base_coord.x() + i][base_coord.y() + j];
				Vector2d grid_v(curr.x(), curr.y());
				double weight = w[i].x() * w[j].y();
				// Velocity
				p.v += weight * grid_v;
				// APIC C, outer product of weighted velocity and dist, paper equation 10
				p.C += 4 * INV_DX * ((weight * grid_v) * dpos.transpose());
			}
		}

		// Advection
		p.x += DT * p.v;

		// MLS-MPM F-update eqn 17
		Matrix2d F = (Matrix2d::Identity() + DT * p.C) * p.F;
		p.F = F;
	}
}

//Rendering
void InitGL(void)
{
	//0.9,0.9,0.9
	//0.5
	glClearColor(0.9f, 0.9f, 0.9f, 1);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(2);
	glMatrixMode(GL_PROJECTION);
}

//Rendering
void Render(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, 0, 1);
	//0.2, 0.6, 1.0 blue
	//0.3, 0.8, 0.3 green
	//1.0, 1.0, 1.0 white
	//glColor4f(0.3f, 0.8f, 0.3f, 1);
	/*
	glBegin(GL_POINTS);
	for (auto& p : particles)
		glVertex2f(p.x(0) * WINDOW_WIDTH, p.x(1) * WINDOW_HEIGHT);
	glEnd();
	*/
	glBegin(GL_POINTS);
	for (auto& p : particles) {
		glColor4f(p.color.x(), p.color.y(), p.color.z(), 1);
		glVertex2f(p.x(0) * WINDOW_WIDTH, p.x(1) * WINDOW_HEIGHT);
	}
	glEnd();
	glColor4f(0.3f, 0.3f, 0.3f, 1);
	for (Rect r : rectangles) {
		glRectf(r.x1 * WINDOW_WIDTH, r.y1 * WINDOW_HEIGHT, r.x2 * WINDOW_WIDTH, r.y2 * WINDOW_HEIGHT);
	}
	glBegin(GL_TRIANGLES);
	for (Tri t : triangles) {
		glVertex2f((float)t.a.x() * WINDOW_WIDTH, (float)t.a.y() * WINDOW_HEIGHT);
		glVertex2f((float)t.b.x() * WINDOW_WIDTH, (float)t.b.y() * WINDOW_HEIGHT);
		glVertex2f((float)t.c.x() * WINDOW_WIDTH, (float)t.c.y() * WINDOW_HEIGHT);
	}
	glEnd();
	
	glutSwapBuffers();
}

void SaveFrame(int frame)
{
	unsigned char* buffer = (unsigned char*)malloc(WINDOW_WIDTH * WINDOW_HEIGHT * 3);
	glReadPixels(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, buffer);
	string filepath = "./imgs/mpm" + to_string(frame) + ".png";
	stbi_flip_vertically_on_write(true);
	stbi_write_png(filepath.c_str(), WINDOW_WIDTH, WINDOW_HEIGHT, 3, buffer, WINDOW_WIDTH * 3);
	frame++;
	free(buffer);
}

void Update(void)
{
	//P2G();
	P2GJello();
	UpdateGridVelocity();
	//cout << "maxparts: " << maxPartsPerGrid << endl;
	//G2P();
	G2PJello();
	step++;
	glutPostRedisplay();
	//Uncomment this to render to file. Causes lag.
	/*
	if (step % INV_FPS == 0) {
		Render();
		SaveFrame(step / INV_FPS);
	}*/
}

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y)
{
	switch (c)
	{
	case ' ':
		if (rectangles.size() > 0) rectangles.clear();
		else addRectangles();
		break;
	case 'r':
	case 'R':
		particles.clear();
		break;
	case 't':
	case 'T':
		if (tetris) tetris = false;
		else tetris = true;
		break;
	case 's':
	case 'S':
		if (triangles.size() > 0) triangles.clear();
		else addSlope();
		break;
	case 'm':
	case 'M':
		InitMPM();
	}
}

void OnClick(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		if (tetris) addTetris(x / (float)WINDOW_WIDTH, 1.f - (y / (float)WINDOW_HEIGHT));
		else addParticles(x / (float)WINDOW_WIDTH, 1.f - (y / (float)WINDOW_HEIGHT));
	}
}

int main(int argc, char** argv)
{
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInit(&argc, argv);
	glutCreateWindow("MPM");
	glutDisplayFunc(Render);
	glutIdleFunc(Update);
	glutKeyboardFunc(Keyboard);
	glutMouseFunc(OnClick);

	InitGL();
	InitMPM();

	glutMainLoop();
	return 0;
}