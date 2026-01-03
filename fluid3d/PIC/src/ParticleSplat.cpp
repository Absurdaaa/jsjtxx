#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../../../third_party/stb/stb_image_write.h"

#include "PIC/include/ParticleSplat.h"
#include <glm/gtc/matrix_transform.hpp>
#include <cmath>
#include <algorithm>

namespace FluidSimulation { namespace PIC3d {

static inline int clampi(int v, int a, int b) { return v < a ? a : (v > b ? b : v); }

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
    ImageRGBA &outImage)
{
    outImage.width = width;
    outImage.height = height;
    outImage.rgba.assign(width * height * 4, 0);

    std::vector<float> density(width * height, 0.0f);

    int radius = std::max(1, (int)std::ceil(3.0f * sigmaPixels));
    int ksize = 2 * radius + 1;
    std::vector<float> kernel(ksize * ksize);
    const float twoSigma2 = 2.0f * sigmaPixels * sigmaPixels;
    float kernelSum = 0.0f;
    for (int dy = -radius; dy <= radius; ++dy) {
        for (int dx = -radius; dx <= radius; ++dx) {
            float r2 = float(dx * dx + dy * dy);
            float w = std::exp(-r2 / twoSigma2);
            kernel[(dy + radius) * ksize + (dx + radius)] = w;
            kernelSum += w;
        }
    }
    if (kernelSum > 0.0f) {
        for (auto &v : kernel) v /= kernelSum;
    }

    for (const ParticleSimple &p : particles) {
        glm::vec4 clip = proj * view * glm::vec4(p.position, 1.0f);
        if (clip.w <= 0.0f) continue;
        glm::vec3 ndc = glm::vec3(clip) / clip.w;
        if (ndc.x < -1.0f - 0.01f || ndc.x > 1.0f + 0.01f || ndc.y < -1.0f - 0.01f || ndc.y > 1.0f + 0.01f) continue;

        float sx = (ndc.x * 0.5f + 0.5f) * float(width);
        float sy = (1.0f - (ndc.y * 0.5f + 0.5f)) * float(height);

        int cx = int(std::floor(sx + 0.5f));
        int cy = int(std::floor(sy + 0.5f));

        int x0 = clampi(cx - radius, 0, width - 1);
        int x1 = clampi(cx + radius, 0, width - 1);
        int y0 = clampi(cy - radius, 0, height - 1);
        int y1 = clampi(cy + radius, 0, height - 1);

        for (int y = y0; y <= y1; ++y) {
            for (int x = x0; x <= x1; ++x) {
                int kx = x - (cx - radius);
                int ky = y - (cy - radius);
                int kidx = ky * ksize + kx;
                float w = kernel[kidx];
                if (w <= 0.0f) continue;
                int idx = y * width + x;
                density[idx] += p.density * w;
            }
        }
    }

    auto toU8 = [](float v)->unsigned char {
        if (v <= 0.0f) return 0;
        if (v >= 1.0f) return 255;
        return (unsigned char)(v * 255.0f + 0.5f);
    };

    for (int i = 0; i < width * height; ++i) {
        float d = density[i];
        float alpha = 1.0f - std::exp(-extinction * d);
        glm::vec3 col = smokeColor * alpha + bgColor * (1.0f - alpha);
        int base = i * 4;
        outImage.rgba[base + 0] = toU8(col.r);
        outImage.rgba[base + 1] = toU8(col.g);
        outImage.rgba[base + 2] = toU8(col.b);
        outImage.rgba[base + 3] = toU8(alpha);
    }

    // write image for quick preview (optional)
    // stbi_write_png requires row stride in bytes
    stbi_write_png("cpu_particle_render.png", width, height, 4, outImage.rgba.data(), width * 4);
}

}}
