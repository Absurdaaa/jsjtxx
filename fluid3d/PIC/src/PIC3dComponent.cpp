/**
 * PIC3dComponent.cpp: 3D PIC 组件实现
 * 实现组件的生命周期管理
 */

#include "PIC3dComponent.h"
#include "Logger.h"
#include "Global.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 关闭组件，释放所有资源
         */
        void PIC3dComponent::shutDown()
        {
            delete renderer;
            delete solver;
            delete grid;
            delete particleSystem;
            renderer = nullptr;
            solver = nullptr;
            grid = nullptr;
            particleSystem = nullptr;
        }

        /**
         * 初始化组件
         * 创建网格、粒子系统、渲染器和求解器
         */
        void PIC3dComponent::init()
        {
            // 如果已存在，先释放
            if (renderer != nullptr || solver != nullptr || grid != nullptr)
            {
                shutDown();
            }

            // 重置计时器
            Glb::Timer::getInstance().clear();

            // 创建网格和粒子系统
            grid = new PICGrid3d();
            particleSystem = new ParticleSystem3d();

            // 记录日志
            Glb::Logger::getInstance().addLog("3D PIC grid created. dimension: " +
                std::to_string(Eulerian3dPara::theDim3d[0]) + "x" +
                std::to_string(Eulerian3dPara::theDim3d[1]) + "x" +
                std::to_string(Eulerian3dPara::theDim3d[2]));

            // 创建渲染器和求解器
            renderer = new Renderer3d();
            renderer->init();
            solver = new Solver3d(*particleSystem, *grid);
        }

        /**
         * 执行一步仿真
         */
        void PIC3dComponent::simulate()
        {
            solver->solve();
        }

        /**
         * 获取渲染结果纹理
         */
        GLuint PIC3dComponent::getRenderedTexture()
        {
            renderer->draw(*particleSystem);
            return renderer->getRenderedTexture();
        }
    }
}
