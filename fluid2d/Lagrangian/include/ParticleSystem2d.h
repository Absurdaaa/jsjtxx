/**
 * ParticleSystem2d.h: 2D粒子系统头文件
 * 定义基于粒子的2D流体仿真数据结构
 */

#pragma once
#ifndef __PARTICAL_SYSTEM_2D_H__
#define __PARTICAL_SYSTEM_2D_H__

#include <vector>
#include <list>
#include <glm/glm.hpp>
#include "Global.h"

#include "Configure.h"

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        /**
         * 2D粒子信息结构体
         * 存储粒子的物理属性
         */
        struct ParticleInfo2d
        {
            alignas(8) glm::vec2 position;      // 位置
            alignas(8) glm::vec2 velocity;      // 速度
            alignas(8) glm::vec2 accleration;    // 加速度
            alignas(4) float density;            // 密度
            alignas(4) float pressure;          // 压力
            alignas(4) float pressDivDens2;      // 压力除以密度的平方
            alignas(4) uint32_t blockId;         // 所属流体块ID
        };

        /**
         * 2D粒子系统类
         * 管理所有粒子的物理属性和空间索引
         */
        class ParticleSystem2d
        {
        public:
            ParticleSystem2d();
            ~ParticleSystem2d();

            // 设置流体容器的大小
            void setContainerSize(glm::vec2 containerCorner, glm::vec2 containerSize);
            int32_t addFluidBlock(glm::vec2 corner, glm::vec2 size, glm::vec2 v0, float particleSpace);
            uint32_t getBlockIdByPosition(glm::vec2 position);
            void updateBlockInfo();

        public:
          // 粒子参数
          float supportRadius = Lagrangian2dPara::supportRadius;       // 支持域半径：影响核函数权重的影响范围（通常用于近邻搜索）
          float supportRadius2 = supportRadius * supportRadius;        // 支持域半径的平方（用于避免在循环中重复计算平方）
          float particleRadius = Lagrangian2dPara::particleRadius;     // 粒子半径（用于碰撞检测或渲染）
          float particleDiameter = Lagrangian2dPara::particleDiameter; // 粒子直径（通常为 2 * particleRadius）
          float particleVolume = 3.1415926f * particleRadius * particleRadius; // 粒子体积（2D中为面积，πr?）

          // 存储全部粒子信息
          std::vector<ParticleInfo2d> particles;

          // 容器参数
          glm::vec2 lowerBound = glm::vec2(FLT_MAX); // 分别表示容器的左下角和右上角坐标
          glm::vec2 upperBound = glm::vec2(-FLT_MAX);
          glm::vec2 containerCenter = glm::vec2(0.0f);
          float scale = Lagrangian2dPara::scale;

          // Block结构（加速临近搜索）
          glm::uvec2 blockNum = glm::uvec2(0);
          glm::vec2 blockSize = glm::vec2(0.0f);
          std::vector<glm::uvec2> blockExtens;
          std::vector<int32_t> blockIdOffs;
        };
    }
}

#endif // !PARTICAL_SYSTEM_H
