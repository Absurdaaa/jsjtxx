/**
 * Configure.h: ϵͳ�����ļ�
 * ����ϵͳ�ĸ������ò�����ȫ�ֱ���
 */

#pragma once
#ifndef __CONFIGURE_H__
#define __CONFIGURE_H__

#include <iostream>
#include <string>
#include "Component.h"
#include <vector>
#include "glm/glm.hpp"

/**
 * ���Բ�ֵ��
 * @param a ��ʼֵ
 * @param b ����ֵ
 * @param t ��ֵ���� [0,1]
 */
#define LERP(a, b, t) (1 - t) * a + t *b

#ifndef __MINMAX_DEFINED
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

// ��Ⱦ�ֱ���
extern int imageWidth;      // ��Ⱦͼ�����
extern int imageHeight;     // ��Ⱦͼ��߶�

// ���ڳߴ�
extern int windowWidth;     // ���ڿ���
extern int windowHeight;    // ���ڸ߶�

// UI����
extern float fontSize;      // �����С

// ����״̬
extern bool simulating;     // �Ƿ����ڷ���

/**
 * 2Dŷ���������������ռ�
 * ����2Dŷ�����������������ò���
 */
namespace Eulerian2dPara
{
    /**
     * ����Դ�ṹ��
     * ��������Դ�ĳ�ʼλ�á��ٶȡ��ܶȺ��¶�
     */
    struct SourceSmoke {
        glm::ivec2 position = glm::ivec2(0);
        glm::vec2 velocity = glm::vec2(0.0f);
        float density = 0.0f;
        float temp = 0.0f;
    };

    extern int theDim2d[];
    extern std::vector<SourceSmoke> source;
    extern float theCellSize2d;
    extern bool addSolid;

    extern float dt;

    extern float contrast;
    extern int drawModel;
    extern int gridNum;


    extern float airDensity;
    extern float ambientTemp;
    extern float boussinesqAlpha;
    extern float boussinesqBeta;
    extern float vorticityConst;
}

/**
 * 3Dŷ���������������ռ�
 * ����3Dŷ�����������������ò���
 */
namespace Eulerian3dPara
{
    /**
     * ����Դ�ṹ��
     * ��������Դ�ĳ�ʼλ�á��ٶȡ��ܶȺ��¶�
     */
    struct SourceSmoke {
        glm::ivec3 position = glm::ivec3(0);
        glm::vec3 velocity = glm::vec3(0.0f);
        float density = 0.0f;
        float temp = 0.0f;
    };

    extern int theDim3d[];
    extern float theCellSize3d;
    extern std::vector<SourceSmoke> source;
    extern bool addSolid;

    extern float contrast;
    extern bool oneSheet;
    extern float distanceX;
    extern float distanceY;
    extern float distanceZ;
    extern bool xySheetsON;
    extern bool yzSheetsON;
    extern bool xzSheetsON;
    extern int xySheetsNum;
    extern int yzSheetsNum;
    extern int xzSheetsNum;
    extern int drawModel;
    extern int gridNumX;
    extern int gridNumY;
    extern int gridNumZ;

    extern float dt;

    extern float airDensity;
    extern float ambientTemp;
    extern float boussinesqAlpha;
    extern float boussinesqBeta;
    extern float vorticityConst;

}

/**
 * 2D PIC 混合方法参数命名空间
 */
namespace PIC2dPara
{
    extern int particlesPerStep;   // 每个源每步发射粒子数
    extern float emissionJitter;   // 发射抖动系数（以cell size为单位）
    extern float wallRestitution;  // 边界弹性系数
}

/**
 * 2D�������շ������������ռ�
 * ����2D�����������������������ò���
 */
namespace Lagrangian2dPara
{
    /**
     * �����ṹ��
     * �����ʼ������λ�á���С����ʼ�ٶȺ����Ӽ��
     */
    struct FluidBlock {
        glm::vec2 lowerCorner = glm::vec2(0.0f, 0.0f);
        glm::vec2 upperCorner = glm::vec2(0.0f, 0.0f);
        glm::vec2 initVel = glm::vec2(0.0f, 0.0f);
        float particleSpace = 0.02f;
    };

    extern float scale;
    extern std::vector<FluidBlock> fluidBlocks;

    extern float dt;
    extern int substep;
    extern float maxVelocity;
    extern float velocityAttenuation;
    extern float eps;

    extern float supportRadius;// �ڴ˰뾶�ڵ��������ӻ�Ե�ǰ���Ӳ���Ӱ��
    extern float particleRadius;// ���Ӱ뾶
    extern float particleDiameter;// ����ֱ��
    extern float gravityX;// x���ϵļ��ٶ�
    extern float gravityY;// y���ϵļ��ٶ�
    extern float density;// �����ܶ�
    extern float stiffness;// �ն�
    extern float exponent;// ѹ�����㹫ʽ�е�ָ��
    extern float viscosity;// ճ��
}

// 2D����������Ȫ��������
namespace Lagrangian2dFountainPara
{
        /**
     * �����ṹ��
     * �����ʼ������λ�á���С����ʼ�ٶȺ����Ӽ��
     */

    extern float scale;

    extern float dt;
    extern int substep;
    extern float maxVelocity;
    extern float velocityAttenuation;
    extern float eps;

    extern float supportRadius;// �ڴ˰뾶�ڵ��������ӻ�Ե�ǰ���Ӳ���Ӱ��
    extern float particleRadius;// ���Ӱ뾶
    extern float particleDiameter;// ����ֱ��
    extern float gravityX;// x���ϵļ��ٶ�
    extern float gravityY;// y���ϵļ��ٶ�
    extern float density;// �����ܶ�
    extern float stiffness;// �ն�
    extern float exponent;// ѹ�����㹫ʽ�е�ָ��
    extern float viscosity;// ճ��
  
  
    extern glm::vec2 containerLower;   // ��Ȫ�������½ǣ���������߽磩
    extern glm::vec2 containerUpper;   // ��Ȫ�������Ͻ�
    
    extern glm::vec2 emitterLower;     // ����������½ǣ������ײ�ʵ��������
    extern glm::vec2 emitterUpper;     // ����������Ͻǣ�������ڿ���/�߶ȣ�
    extern glm::vec2 emitterVelocity;  // ���ӳ�ʼ�ٶȣ�ͨ��ָ�� +Y��
    extern float emitterJitter;        // �ٶȶ���ǿ�ȣ����͹��ȹ����µ����ƣ�
    
    extern float particleSpacing;      // �����ڲ����Ӳ������
    extern int particlesPerStep;       // ÿ֡�������������Ӱ��ˮ���ܶȣ�
    extern size_t maxParticles;        // ����ϵͳ�������ڵ��������
}

/**
 * 3D�������շ������������ռ�
 * ����3D�����������������������ò���
 */
namespace Lagrangian3dPara
{
    /**
     * �����ṹ��
     * �����ʼ������λ�á���С����ʼ�ٶȺ����Ӽ��
     */
    struct FluidBlock {
        glm::vec3 lowerCorner = glm::vec3(0.0f, 0.0f, 0.0f);
        glm::vec3 upperCorner = glm::vec3(0.0f, 0.0f, 0.0f);
        glm::vec3 initVel = glm::vec3(0.0f, 0.0f, 0.0f);
        float particleSpace = 0.02f;
    };

    extern float scale;
    extern std::vector<FluidBlock> fluidBlocks;


    extern float dt;// ʱ�䲽��
    extern int substep;// �Ӳ���
    extern float maxVelocity;// ��������ٶ�
    extern float velocityAttenuation;// ��ײ����ٶ�˥��ϵ��
    extern float eps;// һ����С�ľ��룬���ڱ߽紦��

    extern float supportRadius;// ֧�ְ뾶
    extern float particleRadius;// ���Ӱ뾶
    extern float particleDiameter;// ����ֱ��

    extern float gravityX;// x������
    extern float gravityY;// y������
    extern float gravityZ;// z������

    extern float density;// �����ܶ�
    extern float stiffness;// �ն�
    extern float exponent;//  ѹ�����㹫ʽ�е�ָ��
    extern float viscosity;// ճ��
    
    // XPSH��ز���
    extern float xsph_c; // XSPH����ϵ��
}

// ��Դ·��
extern std::string shaderPath;      // ��ɫ���ļ�·��
extern std::string picturePath;      // ͼƬ�ļ�·��

// �������
extern std::vector<Glb::Component *> methodComponents;  // ���з���������б�

#endif // !__CONFIGURE_H__
