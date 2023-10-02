


struct Particle{
    float m; // mass
    float *x; // position vector 
    float *v; // velocity vector
    float *f; // force accumulator 
};

struct ParticleSystem{
    Particle *p; // array of particles
    int n;      // number of particles
    float t;    // simulation clock 
};

// length of state derivation and force vectors 
int ParticleDims(ParticleSystem* ps){
    return( 6 * ps->n);
}

// gatter state form the particles into dst
int ParticleGetState(ParticleSystem* ps, float *dst){
    
    int i;
    for(i =0;i < ps->n; i++){
        *(dst++) = ps->p[i].x[0]; ///x
        *(dst++) = ps->p[i].x[1]; //y
        *(dst++) = ps->p[i].x[2];//z
        *(dst++) = ps->p[i].v[0]; //vx
        *(dst++) = ps->p[i].v[1]; //vy
        *(dst++) = ps->p[i].v[2]; //vz
    }
}

// scatter state from src into the particles 
int ParticleSetState(ParticleSystem* ps, float *src){
  int i;
  for(i=0; i < ps->n; i++){
    ps->p[i].x[0] = *(src++);
    ps->p[i].x[1] = *(src++);
    ps->p[i].x[2] = *(src++);
    ps->p[i].v[0] = *(src++);
    ps->p[i].v[1] = *(src++);
    ps->p[i].v[2] = *(src++);
  }
}
/* zero the force accumulators */
void Clear_Forces(ParticleSystem* ps){

}
/* magic force function */
void Compute_Forces(ParticleSystem* ps){

}
/* calculate derivative, place in dst */
int ParticleDerivative(ParticleSystem* ps, float *dst){
 int i;
 Clear_Forces(ps);   
 Compute_Forces(ps);
  for(i=0; i < ps->n; i++){
        float m  = ps->p[i].m;
        *(dst++) = ps->p[i].v[0];    /* xdot = v */
        *(dst++) = ps->p[i].v[1];
        *(dst++) = ps->p[i].v[2];
        *(dst++) = ps->p[i].f[0]/m; /* vdot = f/m */
        *(dst++) = ps->p[i].f[1]/m;
        *(dst++) = ps->p[i].f[2]/m;
    } 
}
void ParticleDeriv(ParticleSystem* ps, float * vector){

}
void ScaleVector(float *temp1, float deltaT){

}
void AddVectors(float* a, float *b, float *c){

}

void EulerStep(ParticleSystem* ps, float DeltaT){
    float *temp1,*temp2,*temp3, deltaT;
  ParticleDeriv(ps,temp1);
  ScaleVector(temp1,deltaT);
  ParticleGetState(ps,temp2);
  AddVectors(temp1,temp2,temp2);  /* add -> temp2 */
  ParticleSetState(ps,temp2);      /* update state */
  ps->t += DeltaT;                 /* update time */
}