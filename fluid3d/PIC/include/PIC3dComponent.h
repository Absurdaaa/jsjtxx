/**
 * PIC3dComponent.h: 3D PIC 组件头文件
 * 定义 3D PIC 仿真的组件类，继承自 Component 基类
 */

#pragma once
#ifndef __PIC_3D_COMPONENT_H__
#define __PIC_3D_COMPONENT_H__

#include "Renderer3d.h"
#include "Solver3d.h"
#include "PICGrid3d.h"
#include "ParticleSystem3d.h"

#include "Component.h"
#include "Configure.h"
#include "Global.h"

namespace FluidSimulation
{
    namespace PIC3d
    {
        /**
         * 3D PIC 组件类
         * 管理 PIC 仿真的生命周期：初始化、仿真、渲染、关闭
         */
        class PIC3dComponent : public Glb::Component
        {
        public:
            Renderer3d *renderer;           // 渲染器
            Solver3d *solver;               // 求解器
            PICGrid3d *grid;                // MAC 网格
            ParticleSystem3d *particleSystem;  // 粒子系统

            /**
             * 构造函数
             * @param description 组件描述
             * @param id 组件 ID
             */
            PIC3dComponent(char *description, int id)
            {
                this->description = description;
                this->id = id;
                renderer = nullptr;
                solver = nullptr;
                grid = nullptr;
                particleSystem = nullptr;
            }

            virtual void shutDown();           // 关闭组件，释放资源
            virtual void init();               // 初始化组件
            virtual void simulate();           // 执行一步仿真
            virtual GLuint getRenderedTexture();  // 获取渲染结果
        };
    }
}

#endif
