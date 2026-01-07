#include "Configure.h"

// ϵͳ���ò���
int imageWidth = 600;       // ��Ⱦ�ֱ��ʿ���
int imageHeight = 600;      // ��Ⱦ�ֱ��ʸ߶�

int windowWidth = 1080;     // ����Ĭ�Ͽ���
int windowHeight = 960;     // ����Ĭ�ϸ߶�

float fontSize = 16.0f;     // �����С

bool simulating = false;    // ģ��״̬�������л���ͣ

// 2Dŷ��������ز���
namespace Eulerian2dPara
{
    // MAC�������
    int theDim2d[2] = {100, 100};   // ����ά��
    float theCellSize2d = 0.5;      // ����Ԫ��С
    
    // ����Դ�����ò���
    std::vector<SourceSmoke> source = {
        {   // Դ��λ��                      // ��ʼ�ٶ�           // �ܶ�  // �¶�
            glm::ivec2(theDim2d[0] / 3, 0), glm::vec2(0.0f, 1.0f), 1.0f, 1.0f
        }
    };

    bool addSolid = true;           // �Ƿ����ӹ���߽�

    // ��Ⱦ������
    float contrast = 1;             // �Աȶ�
    int drawModel = 0;             // ����ģʽ
    int gridNum = theDim2d[0];     // ��������

    // ���������
    float dt = 0.01;               // ʱ�䲽��
    float airDensity = 1.3;        // �����ܶ�
    float ambientTemp = 0.0;       // �����¶�
    float boussinesqAlpha = 500.0; // Boussinesq����ʽ�е�alpha����
    float boussinesqBeta = 2500.0; // Boussinesq����ʽ�е�beta����

    // 涡量增强系数（0 表示关闭）
    float vorticityConst = 0.0f;
}

// 3Dŷ��������ز���
namespace Eulerian3dPara
{
    // MAC�������
    // 默认使用正方体场景：x = y = z
    int theDim3d[3] = {36, 36, 36}; // ����ά�� (正方体)
    float theCellSize3d = 0.5;      // ����Ԫ��С
    std::vector<SourceSmoke> source = {
        {glm::ivec3(theDim3d[0] / 2, theDim3d[1] / 2, 0), glm::vec3(0.0f, 0.0f, 1.0f), 1.0f, 1.0f}
    };
    bool addSolid = true;           // �Ƿ����ӹ���߽�

    // ��Ⱦ������
    float contrast = 1;             // �Աȶ�
    bool oneSheet = true;           // �Ƿ�ֻ��ʾһ����Ƭ
    float distanceX = 0.51;         // X����Ƭλ��
    float distanceY = 0.51;         // Y����Ƭλ��
    float distanceZ = 0.985;        // Z����Ƭλ��
    bool xySheetsON = true;         // �Ƿ���ʾXYƽ����Ƭ
    bool yzSheetsON = true;         // �Ƿ���ʾYZƽ����Ƭ
    bool xzSheetsON = true;         // �Ƿ���ʾXZƽ����Ƭ
    int drawModel = 0;              // ����ģʽ
    int gridNumX = (int)((float)theDim3d[0] / theDim3d[2] * 100);  // X��������
    int gridNumY = (int)((float)theDim3d[1] / theDim3d[2] * 100);  // Y��������
    int gridNumZ = 100;             // Z��������
    int xySheetsNum = 3;            // XYƽ����Ƭ����
    int yzSheetsNum = 3;            // YZƽ����Ƭ����
    int xzSheetsNum = 3;            // XZƽ����Ƭ����

    // ���������
    float dt = 0.01;                // ʱ�䲽��
    float airDensity = 1.3;         // �����ܶ�
    float ambientTemp = 0.0;        // �����¶�
    float boussinesqAlpha = 500.0;  // Boussinesq����ʽ�е�alpha����
    float boussinesqBeta = 2500.0;  // Boussinesq����ʽ�е�beta����

    // 涡量增强系数（0 表示关闭）
    float vorticityConst = 0.0f;
}

// 2D PIC 混合方法参数
namespace PIC2dPara
{
  // yuanben
  // MAC�������
  int theDim2d[2] = {100, 100}; // ����ά��
  float theCellSize2d = 0.5;    // ����Ԫ��С

  // ����Դ�����ò���
  // 默认场景：左侧中间喷出，向右流动，以便绕过中间空中障碍
  std::vector<SourceSmoke> source = {
      {// 源位置                      // ��ʼ�ٶ�           // �ܶ�  // �¶�
       glm::ivec2(0, theDim2d[1] / 2), 
       glm::vec2(500.0f, 0.0f), 
       1.0f, 
       1.0f}
       };

  bool addSolid = true; // �Ƿ����ӹ���߽�

  // ��Ⱦ������
  float contrast = 1;        // �Աȶ�
  int drawModel = 0;         // ����ģʽ
  int gridNum = theDim2d[0]; // ��������

  // ���������
  float dt = 0.01;               // ʱ�䲽��
  float airDensity = 1.3;        // �����ܶ�
  float ambientTemp = 0.0;       // �����¶�
  float boussinesqAlpha = 500.0; // Boussinesq����ʽ�е�alpha����
  float boussinesqBeta = 2500.0; // Boussinesq����ʽ�е�beta����

    // 涡量增强系数（0 表示关闭；常用范围 1~50，取决于尺度/步长）
    float vorticityConst = 20.0f;

  int particlesPerStep = 500;
  float emissionJitter = 0.3f;
  float wallRestitution = 0.1f;
    int emitterRadius = 10; // 发射源半径（格子数），1 表示覆盖中心及邻格

    // 风：默认关闭（0）。如果你希望一直向右吹，可以设为正值（例如 50~200）
    float windX = 500.0f;
}

// 3D PIC 混合方法参数
namespace PIC3dPara
{
    // 使用立方体场景，方便可视化与设置
    int theDim3d[3] = {36, 36, 36}; // ����ά�� (正方体)
  float theCellSize3d = 0.5;      // ����Ԫ��С
  std::vector<SourceSmoke> source = {
      {glm::ivec3(0, theDim3d[1] / 2, theDim3d[2] / 2), glm::vec3(100.0f, 0.0f, 0.0f), 1.0f, 1.0f}};
  bool addSolid = true; // �Ƿ����ӹ���߽�

  // ��Ⱦ������
  float contrast = 1;                                           // �Աȶ�
  bool oneSheet = true;                                         // �Ƿ�ֻ��ʾһ����Ƭ
  float distanceX = 0.51;                                       // X����Ƭλ��
  float distanceY = 0.51;                                       // Y����Ƭλ��
  float distanceZ = 0.985;                                      // Z����Ƭλ��
  bool xySheetsON = true;                                       // �Ƿ���ʾXYƽ����Ƭ
  bool yzSheetsON = true;                                       // �Ƿ���ʾYZƽ����Ƭ
  bool xzSheetsON = true;                                       // �Ƿ���ʾXZƽ����Ƭ
  int drawModel = 0;                                            // ����ģʽ
  int gridNumX = (int)((float)theDim3d[0] / theDim3d[2] * 100); // X��������
  int gridNumY = (int)((float)theDim3d[1] / theDim3d[2] * 100); // Y��������
  int gridNumZ = 100;                                           // Z��������
  int xySheetsNum = 3;                                          // XYƽ����Ƭ����
  int yzSheetsNum = 3;                                          // YZƽ����Ƭ����
  int xzSheetsNum = 3;                                          // XZƽ����Ƭ����

  // ���������
  float dt = 0.01;               // ʱ�䲽��
  float airDensity = 1.3;        // �����ܶ�
  float ambientTemp = 0.0;       // �����¶�
  float boussinesqAlpha = 500.0; // Boussinesq����ʽ�е�alpha����
  float boussinesqBeta = 2500.0; // Boussinesq����ʽ�е�beta����

    // 涡量增强系数（0 表示关闭）
    float vorticityConst = 0.0f;
  int pressureIters = 50;        // Default value for pressure iterations

  // 场景：中心球体半径（cell）。0 表示关闭障碍。
  float sphereRadiusCells = 3.5f;

    // 风：默认关闭（0）。正值向 +X 吹。
    float windX = 0.0f;

  int particlesPerStep = 500;
  float emissionJitter = 0.3f;
  float wallRestitution = 0.1f;
  int emitterRadius = 5;
}

// 2D�������շ�����ز���
namespace Lagrangian2dPara
{
    // ����ϵ��
    float scale = 2;

    // ������ʼ����
    std::vector<FluidBlock> fluidBlocks = {
        {   // ���½�����             // ���Ͻ�����           // ��ʼ�ٶ�              // ���Ӽ��
            glm::vec2(-0.4f, -0.4f), glm::vec2(0.4f, 0.4f), glm::vec2(0.0f, 0.0f), 0.02f
        }
    };

    // ���������
    float dt = 0.0016;              // ʱ�䲽��
    int substep = 1;                // �Ӳ���
    float maxVelocity = 10;         // ��������ٶ�
    float velocityAttenuation = 0.7; // ��ײ����ٶ�˥��ϵ��
    float eps = 1e-5;               // һ����С�ľ��룬���ڱ߽紦������ֹ���Ӵ����߽�

    // ����ϵͳ����
    float supportRadius = 0.04;      // �ڴ˰뾶�ڵ��������ӻ�Ե�ǰ���Ӳ���Ӱ��
    float particleRadius = 0.01;     // ���ӵİ뾶
    float particleDiameter = particleRadius * 2.0;  // ����ֱ��
    float gravityX = 0.0f;          // x���ϵļ��ٶ�
    float gravityY = 9.8f;          // y���ϵļ��ٶ�
    
    float density = 1000.0f;        // �����ܶ�
    float stiffness = 70.0f;        // �ն�
    float exponent = 7.0f;          // ѹ�����㹫ʽ�е�ָ��
    float viscosity = 0.03f;        // ճ��
}

namespace Lagrangian2dFountainPara
{
  // ����ϵ��
  float scale = 2;

  // ���������
  float dt = 0.0016;               // ʱ�䲽��
  int substep = 1;                 // �Ӳ���
  float maxVelocity = 10;          // ��������ٶ�
  float velocityAttenuation = 0.0; // ��ײ����ٶ�˥��ϵ��
  float eps = 1e-5;                // һ����С�ľ��룬���ڱ߽紦������ֹ���Ӵ����߽�

  // ����ϵͳ����
  float supportRadius = 0.04;                    // �ڴ˰뾶�ڵ��������ӻ�Ե�ǰ���Ӳ���Ӱ��
  float particleRadius = 0.01;                   // ���ӵİ뾶
  float particleDiameter = particleRadius * 2.0; // ����ֱ��
  float gravityX = 0.0f;                         // x���ϵļ��ٶ�
  float gravityY = 9.8f;                         // y���ϵļ��ٶ�

  float density = 1000.0f; // �����ܶ�
  float stiffness = 40.0f; // �ն�
  float exponent = 5.0f;   // ѹ�����㹫ʽ�е�ָ��
  float viscosity = 1.0f; // ճ��

  glm::vec2 containerLower = glm::vec2(-1.0f, -1.0f);
  glm::vec2 containerUpper = glm::vec2(1.0f, 1.0f);

  // ��ڸĵ��ױ����룬��ȸ���
  glm::vec2 emitterLower = glm::vec2(-0.05f, -0.98f);
  glm::vec2 emitterUpper = glm::vec2(0.05f, -0.94f);
  
  glm::vec2 emitterVelocity = glm::vec2(0.0f, 8.0f);
  float emitterJitter = 0.4f;

  float particleSpacing = 0.02f;
  int particlesPerStep = 5;
  size_t maxParticles = 30000;
}

// 3D�������շ�����ز���
namespace Lagrangian3dPara
{
    // ����ϵ��
    float scale = 1.2;
    
    // ������ʼ����
    std::vector<FluidBlock> fluidBlocks = {
        {
            glm::vec3(0.05, 0.05, 0.3), glm::vec3(0.45, 0.45, 0.7), glm::vec3(0.0, 0.0, -1.0), 0.02f
        },
        {
            glm::vec3(0.45, 0.45, 0.3), glm::vec3(0.85, 0.85, 0.7), glm::vec3(0.0, 0.0, -1.0), 0.02f
        }   
    };
    
    // ���������
    float dt = 0.002;               // ʱ�䲽��
    int substep = 1;                // �Ӳ���
    float maxVelocity = 10;         // ��������ٶ�
    float velocityAttenuation = 0.7; // ��ײ����ٶ�˥��ϵ��
    float eps = 1e-5;               // һ����С�ľ��룬���ڱ߽紦��

    // ����ϵͳ����
    float supportRadius = 0.04;      // ֧�ְ뾶
    float particleRadius = 0.01;     // ���Ӱ뾶
    float particleDiameter = particleRadius * 2.0;  // ����ֱ��

    float gravityX = 0.0f;          // x������
    float gravityY = 0.0f;          // y������
    float gravityZ = 9.8f;          // z������

    float density = 1000.0f;        // �ܶ�
    float stiffness = 20.0f;        // �ն�
    float exponent = 7.0f;          // ѹ��ָ��
    float viscosity = 8e-5f;        // ճ��
    
    float xsph_c = 0.2f;            // XSPH����ϵ��
}

// �洢ϵͳ���������
std::vector<Glb::Component *> methodComponents;

// ��Դ·��
std::string shaderPath = "../../../../code/resources/shaders";
std::string picturePath = "../../../../code/resources/pictures";