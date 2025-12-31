/**
 * Lagrangian2dFountainComponent.cpp: 2D喷泉场景组件
 */

#include "Lagrangian2dFountainComponent.h"

#include <algorithm>
#include <vector>

namespace FluidSimulation
{
    namespace Lagrangian2d
    {
        void Lagrangian2dFountainComponent::shutDown()
        {
            delete renderer, solver, ps;
            renderer = NULL;
            solver = NULL;
            ps = NULL;
        }

        void Lagrangian2dFountainComponent::init()
        {
            if (renderer != NULL || solver != NULL || ps != NULL)
            {
                shutDown();
            }

            Glb::Timer::getInstance().clear();

            renderer = new Renderer();
            renderer->init();

            ps = new ParticleSystem2d();
            ps->setContainerSize(Lagrangian2dFountainPara::containerLower, Lagrangian2dFountainPara::containerUpper);
            ps->updateBlockInfo();

            Glb::Logger::getInstance().addLog("2d Fountain particle system ready.");

            solver = new Solver(*ps);
            solver->velocityAttenuation_ = Lagrangian2dFountainPara::velocityAttenuation;
            solver->viscosity_ = Lagrangian2dFountainPara::viscosity;
            solver->substeps_ = Lagrangian2dFountainPara::substep;
            solver->dt_ = Lagrangian2dFountainPara::dt;
            solver->density_ = Lagrangian2dFountainPara::density;
            solver->stiffness_ = Lagrangian2dFountainPara::stiffness;
            solver->exponent_ = Lagrangian2dFountainPara::exponent;
            solver->gravityX_ = Lagrangian2dFountainPara::gravityX;
            solver->gravityY_ = Lagrangian2dFountainPara::gravityY;
            solver->maxVelocity_ = Lagrangian2dFountainPara::maxVelocity;
        }

        void Lagrangian2dFountainComponent::simulate()
        {
            emitParticles();

            const int substeps = (Lagrangian2dFountainPara::substep > 0) ? Lagrangian2dFountainPara::substep : 1;
            for (int i = 0; i < substeps; ++i)
            {
                ps->updateBlockInfo();
                solver->solve();
            }

            removeSettledParticles();
        }

        GLuint Lagrangian2dFountainComponent::getRenderedTexture()
        {
            renderer->draw(*ps);
            return renderer->getRenderedTexture();
        }

        void Lagrangian2dFountainComponent::emitParticles()
        {
            if (ps == NULL)
            {
                return;
            }

            const size_t maxParticles = Lagrangian2dFountainPara::maxParticles;
            if (ps->particles.size() >= maxParticles)
            {
                return;
            }

            int spawnCount = Lagrangian2dFountainPara::particlesPerStep > 0 ? Lagrangian2dFountainPara::particlesPerStep : 0;
            spawnCount = std::min<int>(spawnCount, static_cast<int>(maxParticles - ps->particles.size()));

            std::vector<ParticleInfo2d> newParticles;
            newParticles.reserve(spawnCount);

            for (int i = 0; i < spawnCount; ++i)
            {
                ParticleInfo2d particle = buildParticle();
                if (particle.blockId != UINT32_MAX)
                {
                    newParticles.push_back(particle);
                }
            }

            ps->particles.insert(ps->particles.end(), newParticles.begin(), newParticles.end());
        }

        ParticleInfo2d Lagrangian2dFountainComponent::buildParticle()
        {
            ParticleInfo2d p{};

            glm::vec2 emitterMin = Lagrangian2dFountainPara::emitterLower * Lagrangian2dFountainPara::scale;
            glm::vec2 emitterMax = Lagrangian2dFountainPara::emitterUpper * Lagrangian2dFountainPara::scale;

            float rx = mRand.GetUniformRandom();
            float ry = mRand.GetUniformRandom();

            p.position = glm::vec2(
                emitterMin.x + rx * (emitterMax.x - emitterMin.x),
                emitterMin.y + ry * (emitterMax.y - emitterMin.y));

            glm::vec2 jitter(
                (mRand.GetUniformRandom() - 0.5f) * Lagrangian2dFountainPara::emitterJitter,
                (mRand.GetUniformRandom() - 0.5f) * Lagrangian2dFountainPara::emitterJitter);

            p.velocity = Lagrangian2dFountainPara::emitterVelocity + jitter;
            p.accleration = glm::vec2(0.0f);
            p.density = Lagrangian2dFountainPara::density;
            p.pressure = 0.0f;
            p.pressDivDens2 = 0.0f;
            p.blockId = ps->getBlockIdByPosition(p.position);

            return p;
        }

        void Lagrangian2dFountainComponent::removeSettledParticles()
        {
            if (ps == NULL)
            {
                return;
            }

            const float floorY = Lagrangian2dFountainPara::containerLower.y * Lagrangian2dFountainPara::scale + Lagrangian2dFountainPara::eps;
            auto &particles = ps->particles;
            particles.erase(
                std::remove_if(
                    particles.begin(),
                    particles.end(),
                    [&](const ParticleInfo2d &p)
                    {
                        return p.position.y <= floorY && p.velocity.y <= 0.0f;
                    }),
                particles.end());
        }
    }
}
