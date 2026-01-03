#pragma once
#ifndef __PARTICLE_SPLAT_H__
#define __PARTICLE_SPLAT_H__

#include <vector>
#include <glm/glm.hpp>

namespace FluidSimulation { namespace PIC3d {

struct ParticleSimple {
    glm::vec3 position;
    float density;
};

struct ImageRGBA {
    int width = 0;
    int height = 0;
    std::vector<unsigned char> rgba; // width*height*4
};

// Render particles (position,density) to RGBA image using Gaussian splatting on CPU.
void renderParticlesSplatCPU(
    const std::vector<ParticleSimple>& particles,
    const glm::mat4& view,
    const glm::mat4& proj,
    int width,
    int height,
    float sigmaPixels,
    float extinction,
    const glm::vec3& smokeColor,
    const glm::vec3& bgColor,
    ImageRGBA &outImage);

}}

#endif
