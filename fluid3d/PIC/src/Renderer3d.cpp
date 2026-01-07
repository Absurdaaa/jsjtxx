/**
 * Renderer3d.cpp: 3D PIC 渲染器实现
 * 目标：仅保留 Volume 体渲染：将密度上传为 3D texture，并在片元着色器中做 ray marching。
 */

#include "PIC/include/Renderer3d.h"

// Configure.h 中定义的 min/max 宏会污染 <algorithm>（破坏 std::max/std::min/std::clamp）。
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

#include "Camera.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        Renderer3d::Renderer3d() {}

        Renderer3d::~Renderer3d()
        {
            if (RBO) glDeleteRenderbuffers(1, &RBO);
            if (textureID) glDeleteTextures(1, &textureID);
            if (FBO) glDeleteFramebuffers(1, &FBO);

            if (densityTex3D) glDeleteTextures(1, &densityTex3D);
            if (volumeVAO) glDeleteVertexArrays(1, &volumeVAO);

            delete volumeShader;
            volumeShader = nullptr;
            delete container;
            container = nullptr;
        }

        void Renderer3d::init()
        {
            initGLResources();
        }

        void Renderer3d::initGLResources()
        {
            if (inited)
                return;

            // 尺寸：沿用 Eulerian3d 的归一化容器，保证默认相机能看到。
            const int nz = (Eulerian3dPara::theDim3d[2] <= 0) ? 1 : Eulerian3dPara::theDim3d[2];
            const float xSize = (float)Eulerian3dPara::theDim3d[0] / (float)nz;
            const float ySize = (float)Eulerian3dPara::theDim3d[1] / (float)nz;

            container = new Glb::Container();
            container->resetSize(xSize, ySize, 1.0f);
            container->init();

            // volume ray-march shader (true 3D). drawModel==2 时使用。
            volumeShader = new Glb::Shader();
            std::string v2 = shaderPath + "/VolumeRaymarch3d.vert";
            std::string f2 = shaderPath + "/VolumeRaymarch3d.frag";
            volumeShader->buildFromFile(v2, f2);

            // fullscreen triangle (no VBO needed)
            glGenVertexArrays(1, &volumeVAO);

            // FBO 输出
            glGenFramebuffers(1, &FBO);
            glBindFramebuffer(GL_FRAMEBUFFER, FBO);

            glGenTextures(1, &textureID);
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, imageWidth, imageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
            glBindTexture(GL_TEXTURE_2D, 0);

            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureID, 0);

            glGenRenderbuffers(1, &RBO);
            glBindRenderbuffer(GL_RENDERBUFFER, RBO);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, imageWidth, imageHeight);
            glBindRenderbuffer(GL_RENDERBUFFER, 0);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);

            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            glViewport(0, 0, imageWidth, imageHeight);

            inited = true;
        }

        void Renderer3d::ensureDensityTexture3D(PICGrid3d &grid)
        {
            const int nx = grid.dim[PICGrid3d::X];
            const int ny = grid.dim[PICGrid3d::Y];
            const int nz = grid.dim[PICGrid3d::Z];

            if (nx <= 0 || ny <= 0 || nz <= 0)
                return;

            const bool needAlloc = (!densityTex3D) || (nx != densityNx) || (ny != densityNy) || (nz != densityNz);
            if (!needAlloc)
                return;

            densityNx = nx;
            densityNy = ny;
            densityNz = nz;
            densityData.assign((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);

            if (!densityTex3D)
                glGenTextures(1, &densityTex3D);
            glBindTexture(GL_TEXTURE_3D, densityTex3D);
            glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, densityNx, densityNy, densityNz, 0, GL_RED, GL_FLOAT, nullptr);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
            glBindTexture(GL_TEXTURE_3D, 0);
        }

        void Renderer3d::uploadDensityTexture3D(PICGrid3d &grid)
        {
            ensureDensityTexture3D(grid);
            if (!densityTex3D)
                return;

            const int nx = densityNx;
            const int ny = densityNy;
            const int nz = densityNz;

            // pack density in x-major order: idx = (k*ny + j)*nx + i
            for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i)
                    {
                        const float d = (float)grid.mD(i, j, k);
                        densityData[(size_t)((k * ny + j) * nx + i)] = (d < 0.0f) ? 0.0f : d;
                    }

            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
            glBindTexture(GL_TEXTURE_3D, densityTex3D);
            glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, nx, ny, nz, GL_RED, GL_FLOAT, densityData.data());
            glBindTexture(GL_TEXTURE_3D, 0);
        }

        void Renderer3d::draw(const ParticleSystem3d & /*ps*/, PICGrid3d &grid)
        {
            if (!inited)
                initGLResources();

            glBindFramebuffer(GL_FRAMEBUFFER, FBO);
            glViewport(0, 0, imageWidth, imageHeight);
            glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glDisable(GL_DEPTH_TEST);
            glDepthMask(GL_FALSE);

            // true 3D volume rendering (ray-march on 3D density texture)
            uploadDensityTexture3D(grid);

            glDisable(GL_BLEND);

            const int nzNorm = (grid.dim[PICGrid3d::Z] <= 0) ? 1 : grid.dim[PICGrid3d::Z];
            const float xSize = (float)grid.dim[PICGrid3d::X] / (float)nzNorm;
            const float ySize = (float)grid.dim[PICGrid3d::Y] / (float)nzNorm;

            glm::mat4 view = Glb::Camera::getInstance().GetView();
            glm::mat4 projection = Glb::Camera::getInstance().GetProjection();
            glm::mat4 invViewProj = glm::inverse(projection * view);
            glm::vec3 camPos = Glb::Camera::getInstance().GetPosition();

            // marching step: roughly one cell in normalized space (z extent normalized to 1)
            const int maxDim = (densityNx > densityNy) ? ((densityNx > densityNz) ? densityNx : densityNz)
                                                     : ((densityNy > densityNz) ? densityNy : densityNz);
            const float stepSize = (maxDim > 0) ? (1.0f / (float)maxDim) : (1.0f / 64.0f);

            volumeShader->use();
            volumeShader->setMat4("invViewProj", invViewProj);
            volumeShader->setVec3("cameraPos", camPos);
            volumeShader->setVec3("boxMin", glm::vec3(0.0f, 0.0f, 0.0f));
            volumeShader->setVec3("boxMax", glm::vec3(xSize, ySize, 1.0f));
            volumeShader->setVec2("viewport", (float)imageWidth, (float)imageHeight);
            volumeShader->setFloat("stepSize", stepSize);
            // 用 contrast 同时作为“密度强度/亮度”调节入口
            volumeShader->setFloat("densityScale", PIC3dPara::contrast);
            volumeShader->setVec3("background", glm::vec3(0.05f, 0.05f, 0.05f));
            volumeShader->setVec3("smokeColor", glm::vec3(1.0f, 1.0f, 1.0f));

            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_3D, densityTex3D);
            glUniform1i(glGetUniformLocation(volumeShader->getId(), "densityTex"), 0);

            glBindVertexArray(volumeVAO);
            glDrawArrays(GL_TRIANGLES, 0, 3);
            glBindVertexArray(0);

            if (container)
                container->draw();

            glDepthMask(GL_TRUE);
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

        GLuint Renderer3d::getRenderedTexture()
        {
            return textureID;
        }
    }
}
