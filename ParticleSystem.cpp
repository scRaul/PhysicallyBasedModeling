#include <iostream>
#include <vector>

class Vector3
{
public:
    float x, y, z;

    Vector3(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}

    // Addition operator
    Vector3 operator+(const Vector3 &other) const
    {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }
      Vector3 operator+=(const Vector3 &other) const
    {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    // Subtraction operator
    Vector3 operator-(const Vector3 &other) const
    {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    // Scalar multiplication
    Vector3 operator*(float scalar) const
    {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }
};

class Particle
{
public:
    Vector3 position;
    Vector3 velocity;
    Vector3 forceAccumulator;
    float mass;
    float time;

    Particle(float mass) : mass(mass), time(0.0f)
    {
        position = Vector3(0.0f, 0.0f, 0.0f);
        velocity = Vector3(0.0f, 0.0f, 0.0f);
        forceAccumulator = Vector3(0.0f, 0.0f, 0.0f);
    }
};

class Wall
{
public:
    Vector3 position;
    Vector3 normal;

    Wall(const Vector3 &position, const Vector3 &normal) : position(position), normal(normal) {}
};

class IForce
{
public:
    virtual void apply_force() = 0;
};

class GravityForce : public IForce
{
private:
    Particle *particle;
    float gravity;

public:
    GravityForce(Particle *particle, float gravity) : particle(particle), gravity(gravity) {}

    void apply_force() override
    {
        particle->forceAccumulator += Vector3(0.0f, -particle->mass * gravity, 0.0f);
    }
};

class ParticleSystem
{

public:
    std::vector<Particle *> particles;
    std::vector<IForce *> forces;
    std::vector<Wall> walls;
    float time;
    ParticleSystem() : time(0.0f) {}

    void add_particle(Particle *particle)
    {
        particles.push_back(particle);
    }

    void add_force(IForce *force)
    {
        forces.push_back(force);
    }

    void add_wall(const Vector3 &position, const Vector3 &normal)
    {
        walls.emplace_back(position, normal);
    }
};
// length of state derivation and force vectors
int ParticleDims(ParticleSystem *ps)
{
    return (6 * ps->particles.size());
}

// gatter state form the particles into dst
int ParticleGetState(ParticleSystem *ps, float *dst)
{

    int i;
    for (i = 0; i < ps->particles.size(); i++)
    {
        *(dst++) = ps->particles[i]->position.x;
        *(dst++) = ps->particles[i]->position.y;
        *(dst++) = ps->particles[i]->position.z;
        *(dst++) = ps->particles[i]->velocity.x;
        *(dst++) = ps->particles[i]->velocity.x;
        *(dst++) = ps->particles[i]->velocity.x;
    }
}

// scatter state from src into the particles
int ParticleSetState(ParticleSystem *ps, float *src)
{
    int i;
    for (i = 0; i < ps->particles.size(); i++)
    {
        ps->particles[i]->position.x = *(src++);
        ps->particles[i]->position.y = *(src++);
        ps->particles[i]->position.z = *(src++);
        ps->particles[i]->velocity.x = *(src++);
        ps->particles[i]->velocity.y = *(src++);
        ps->particles[i]->velocity.z = *(src++);
    }
}
void Clear_Forces(ParticleSystem *ps)
{
    int i = 0;
    for (i = 0; i < ps->particles.size(); i++)
    {
        ps->particles[i]->forceAccumulator = 0;
    }
}
void Compute_Forces(ParticleSystem *ps)
{
    int i = 0;
    for (i = 0; i < ps->forces.size(); i++)
    {
        ps->forces[i]->apply_force();
    }
}
/* calculate derivative, place in dst */
int ParticleDerivative(ParticleSystem *ps, float *dst)
{
    int i;
    Clear_Forces(ps);
    Compute_Forces(ps);
    for (i = 0; i < ps->particles.size(); i++)
    {
        float m = ps->particles[i]->mass;
        *(dst++) = ps->particles[i]->velocity.x;
        *(dst++) = ps->particles[i]->velocity.y;
        *(dst++) = ps->particles[i]->velocity.z;
        *(dst++) = ps->particles[i]->forceAccumulator.x / m;
        *(dst++) = ps->particles[i]->forceAccumulator.y / m;
        *(dst++) = ps->particles[i]->forceAccumulator.z / m;
    }
}

void ScaleVector(float* arr,int size,float scale){
    for(int i = 0; i < size; i++){
        arr[i] *= scale;
    }
}

void AddVectors(int size, float *a, float* b, float *c){
    for(int i = 0; i < size; i++){
        c[i] = a[i] + b[i];
    }
}

void EulerStep(ParticleSystem *ps, float deltaT)
{
    int dim = ParticleDims(ps);
    float *temp1 = new float[dim];
    float *temp2 = new float[dim];
    ParticleDerivative(ps, temp1);
    ScaleVector(temp1,dim,deltaT);
    ParticleGetState(ps, temp2);
    AddVectors(dim,temp1,temp2,temp2);  /* add -> temp2 */
    ParticleSetState(ps, temp2); /* update state */
    ps->time += deltaT;          /* update time */
}
int main()
{
    Particle particle(1.0f); // Mass of 1.0
    particle.position = Vector3(0, 0, 0);
    GravityForce gravityForce(&particle, 9.81f); // Standard gravity

    ParticleSystem particleSystem;
    particleSystem.add_particle(&particle);
    particleSystem.add_force(&gravityForce);

    // Add walls
    Vector3 wallPosition(0.0f, 0.0f, 0.0f);
    Vector3 wallNormal(0.0f, 1.0f, 0.0f); // Wall is on the ground plane
    particleSystem.add_wall(wallPosition, wallNormal);

    float deltaT = 0.1f;

    for(float i = 0.0f; i < 1.1f; i+= deltaT){
        EulerStep(&particleSystem,deltaT);
        std::cout << "Time: " << particleSystem.time << "s, Position: (" << particle.position.x << ", " << particle.position.y << ", " << particle.position.z << ")\n";
    }

    return 0;
}
