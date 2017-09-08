
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// core/scene.cpp*
#include "scene.h"
#include "camera.h"
#include "film.h"
#include "sampler.h"
#include "parallel.h"
#include "progressreporter.h"
#include "stats.h"
#include <boost/numeric/odeint.hpp>

namespace pbrt {

STAT_COUNTER("Intersections/Regular ray intersection tests",
             nIntersectionTests);
STAT_COUNTER("Intersections/Shadow ray intersection tests", nShadowTests);

// Scene Method Definitions
bool Scene::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ++nIntersectionTests;
    DCHECK_NE(ray.d, Vector3f(0,0,0));
    return aggregate->Intersect(ray, isect);
    //return IntersectCu(ray, isect);
}

bool Scene::IntersectP(const Ray &ray) const {
    ++nShadowTests;
    DCHECK_NE(ray.d, Vector3f(0,0,0));
    return aggregate->IntersectP(ray);
}

bool Scene::IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *isect,
                        Spectrum *Tr) const {
    *Tr = Spectrum(1.f);
    while (true) {
        bool hitSurface = Intersect(ray, isect);
        // Accumulate beam transmittance for ray segment
        if (ray.medium) *Tr *= ray.medium->Tr(ray, sampler);

        // Initialize next ray segment or terminate transmittance computation
        if (!hitSurface) return false;
        if (isect->primitive->GetMaterial() != nullptr) return true;
        ray = isect->SpawnRay(ray.d);
    }
}

int Scene::forward(const Ray &ray, state_type &prev, state_type &state, SurfaceInteraction *isect) const{
   Vector3f vc1 = Vector3f(prev[0], prev[1], prev[2]);
   Vector3f vc2 = Vector3f(state[0],state[1],state[2]);

   return aggregate->IntersectCu(ray, vc1, vc2, isect);
}

bool Scene::IntersectCu(const Ray &ray, SurfaceInteraction *isect) const {
    ++nIntersectionTests;
    DCHECK_NE(ray.d, Vector3f(0,0,0));

    state_type x(6);
    x[0] = ray.o.x;
    x[1] = ray.o.y;
    x[2] = ray.o.z;

    x[3] = ray.d.x;
    x[4] = ray.d.y;
    x[5] = ray.d.z;

    boost::numeric::odeint::runge_kutta4< state_type > rk;

    double dt = DT_INIT;
    double t = 0.0;

    state_type prev = x;
    for(int i = 0; i<ITERATIONS; ++i, t+=dt ){
        prev[0] = x[0];
        prev[1] = x[1];
        prev[2] = x[2];
        prev[3] = x[3];
        prev[4] = x[4];
        prev[5] = x[5];
        rk.do_step(ray, x, t, dt);
        int stage = forward(ray, prev, x, isect);
        if(stage==1)
            break;
        else if(stage == 0)
        {
            x[0] = prev[0];
            x[1] = prev[1];
            x[2] = prev[2];
            x[3] = prev[3];
            x[4] = prev[4];
            x[5] = prev[5];
            dt = dt/2.0;
            t -= dt;
        }
    }
    if(isect->primitive==nullptr)
      return false;
    else{
      /*
      o.x = x[0];
      o.y = x[1];
      o.z = x[2];

      d.x = x[3];
      d.y = x[4];
      d.z = x[5];
      d = d.norm();
      */
      return true;
    }
}

}  // namespace pbrt
