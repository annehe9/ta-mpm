#if __APPLE__
#include <GLUT/glut.h>
#else
#include <windows.h>
#include <GL/glut.h>
#endif

#include <iostream>
#include <vector>
using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// "Particle-Based Fluid Simulation for Interactive Applications"
// solver parameters
const static Vector2f G(0.f, 12000 * -9.8f); // external (gravitational) forces
const static float REST_DENS = 1000.f;		 // rest density
const static float GAS_CONST = 2000.f;		 // const for equation of state
const static float H = 16.f;				 // kernel radius
const static float HSQ = H * H;				 // radius^2 for optimization
const static float MASS = 65.f;				 // assume all particles have the same mass
const static float VISC = 250.f;			 // viscosity constant
const static float DT = 0.0001f;			 // integration timestep

// smoothing kernels defined in Müller and their gradients
const static float POLY6 = 315.f / (65.f * M_PI * pow(H, 9.f));
const static float SPIKY_GRAD = -45.f / (M_PI * pow(H, 6.f));
const static float VISC_LAP = 45.f / (M_PI * pow(H, 6.f));

// simulation parameters
const static float EPS = H; // boundary epsilon
const static float BOUND_DAMPING = -0.5f;

// particle data structure
// stores position, velocity, and force for integration
// stores density (rho) and pressure values for SPH
struct Particle
{
	Particle(float _x, float _y, float m) : x(_x, _y), v(0.f, 0.f), mass(m) {}
	Vector2f x, v;
	float mass;
};

struct Cell
{
	Cell() : v(0.f, 0.f), mass(0.f) {}
	Vector2f v;
	float mass;

};

// solver data
static vector<Particle> particles;
static vector<Cell> grid;

// interaction
const static int MAX_PARTICLES = 2500;
const static int DAM_PARTICLES = 500;
const static int BLOCK_PARTICLES = 250;
const static int NUM_PARTICLES = 500;
const static int GRID_RES = 64;
const static int NUM_CELLS = GRID_RES * GRID_RES;

// rendering projection parameters
const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
const static double VIEW_WIDTH = 1.5 * 800.f;
const static double VIEW_HEIGHT = 1.5 * 600.f;

void InitMPM(void)
{
	cout << "initializing mpm with " << DAM_PARTICLES << " particles" << endl;
	for (float y = EPS; y < VIEW_HEIGHT - EPS * 2.f; y += H) {
		for (float x = VIEW_WIDTH / 4; x <= VIEW_WIDTH / 2; x += H) {
			if (particles.size() < DAM_PARTICLES)
			{
				float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
				particles.push_back(Particle(x + jitter, y, 1.f));
			}
		}
	}
		
	for (int i = 0; i < NUM_CELLS; i++) {
		grid.push_back(Cell());
	}
}

void Integrate(void)
{
	for (auto& p : particles)
	{
		// forward Euler integration
		p.v += DT * p.f / p.rho;
		p.x += DT * p.v;

		// enforce boundary conditions
		if (p.x(0) - EPS < 0.0f)
		{
			p.v(0) *= BOUND_DAMPING;
			p.x(0) = EPS;
		}
		if (p.x(0) + EPS > VIEW_WIDTH)
		{
			p.v(0) *= BOUND_DAMPING;
			p.x(0) = VIEW_WIDTH - EPS;
		}
		if (p.x(1) - EPS < 0.0f)
		{
			p.v(1) *= BOUND_DAMPING;
			p.x(1) = EPS;
		}
		if (p.x(1) + EPS > VIEW_HEIGHT)
		{
			p.v(1) *= BOUND_DAMPING;
			p.x(1) = VIEW_HEIGHT - EPS;
		}
	}
}

void ComputeDensityPressure(void)
{
	for (auto& pi : particles)
	{
		pi.rho = 0.f;
		for (auto& pj : particles)
		{
			Vector2f rij = pj.x - pi.x;
			float r2 = rij.squaredNorm();

			if (r2 < HSQ)
			{
				// this computation is symmetric
				pi.rho += MASS * POLY6 * pow(HSQ - r2, 3.f);
			}
		}
		pi.p = GAS_CONST * (pi.rho - REST_DENS);
	}
}

void ComputeForces(void)
{
	for (auto& pi : particles)
	{
		Vector2f fpress(0.f, 0.f);
		Vector2f fvisc(0.f, 0.f);
		for (auto& pj : particles)
		{
			if (&pi == &pj)
				continue;

			Vector2f rij = pj.x - pi.x;
			float r = rij.norm();

			if (r < H)
			{
				// compute pressure force contribution
				fpress += -rij.normalized() * MASS * (pi.p + pj.p) / (2.f * pj.rho) * SPIKY_GRAD * pow(H - r, 2.f);
				// compute viscosity force contribution
				fvisc += VISC * MASS * (pj.v - pi.v) / pj.rho * VISC_LAP * (H - r);
			}
		}
		Vector2f fgrav = G * pi.rho;
		pi.f = fpress + fvisc + fgrav;
	}
}

void ClearGrid(void)
{
	for (Cell c : grid) {
		c.mass = 0;
		c.v = Vector2f::Zero();
	}
}

void P2G(void)
{
	for (Particle p : particles) {
		int cell_idx = p.x(0) + (p.x(1) * VIEW_WIDTH);

	}
}

void Update(void)
{
	ComputeDensityPressure();
	ComputeForces();
	Integrate();

	glutPostRedisplay();
}


//Rendering
void InitGL(void)
{
	glClearColor(0.9f, 0.9f, 0.9f, 1);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(H / 2.f);
	glMatrixMode(GL_PROJECTION);
}

//Rendering
void Render(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	glLoadIdentity();
	glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);

	glColor4f(0.2f, 0.6f, 1.0f, 1);
	glBegin(GL_POINTS);
	for (auto& p : particles)
		glVertex2f(p.x(0), p.x(1));
	glEnd();

	glutSwapBuffers();
}

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y)
{
	switch (c)
	{
	case ' ':
		if (particles.size() >= MAX_PARTICLES)
			std::cout << "maximum number of particles reached" << std::endl;
		else
		{
			unsigned int placed = 0;
			for (float y = VIEW_HEIGHT / 1.5f - VIEW_HEIGHT / 5.f; y < VIEW_HEIGHT / 1.5f + VIEW_HEIGHT / 5.f; y += H * 0.95f)
				for (float x = VIEW_WIDTH / 2.f - VIEW_HEIGHT / 5.f; x <= VIEW_WIDTH / 2.f + VIEW_HEIGHT / 5.f; x += H * 0.95f)
					if (placed++ < BLOCK_PARTICLES && particles.size() < MAX_PARTICLES)
						particles.push_back(Particle(x, y));
		}
		break;
	case 'r':
	case 'R':
		particles.clear();
		InitMPM();
		break;
	}
}

int main(int argc, char** argv)
{
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInit(&argc, argv);
	glutCreateWindow("Muller SPH");
	glutDisplayFunc(Render);
	glutIdleFunc(Update);
	glutKeyboardFunc(Keyboard);

	InitGL();
	InitMPM();

	glutMainLoop();
	return 0;
}