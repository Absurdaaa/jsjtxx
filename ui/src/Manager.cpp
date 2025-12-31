/**
 * Manager.cpp: UI������ʵ���ļ�
 * ʵ��UIϵͳ�Ĺ�������ͼ��ʾ
 */

#include "Manager.h"

namespace FluidSimulation
{
    /**
     * ��ʼ��������
     * ����������ͼ��ע������ģ�����
     * @param window GLFW����
     */
    void Manager::init(GLFWwindow* window) {
        // ���洰��ָ��
        this->window = window;

        // ����������ͼ
        inspectorView = new InspectorView(window);
        projectView =  new ProjectView(window);
        sceneView = new SceneView(window);

        // ������ע����������ģ�����
        int id = 0;
        methodComponents.push_back(new Lagrangian2d::Lagrangian2dComponent("Lagrangian 2d", id++));
        methodComponents.push_back(new Lagrangian2d::Lagrangian2dFountainComponent("Lagrangian 2d Fountain", id++));
        methodComponents.push_back(new Eulerian2d::Eulerian2dComponent("Eulerian 2d", id++));
        methodComponents.push_back(new PIC2d::PICComponent("PIC 2d", id++));
        methodComponents.push_back(new Lagrangian3d::Lagrangian3dComponent("Lagrangian 3d", id++));
        methodComponents.push_back(new Eulerian3d::Eulerian3dComponent("Eulerian 3d", id++));
        methodComponents.push_back(new PIC3d::PIC3dComponent("PIC 3d", id++));
        // TODO(optional): ���Ӹ���ģ�ⷽ��
    }

    /**
     * ��ʾ������ͼ
     * ����ͣ���ռ䲢��ʾ������ͼ����������ͼ����Ŀ��ͼ
     */
	void Manager::displayViews() {
        // ����ͣ���ռ��־
        ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;
        ImGui::DockSpaceOverViewport(nullptr, dockspace_flags);

        // ����������ɫ
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.5f, 0.5f, 0.5f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_BorderShadow, ImVec4(0.5f, 0.5f, 0.5f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.1f, 0.1f, 0.1f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_TitleBg, ImVec4(0.2f, 0.2f, 0.2f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_TitleBgActive, ImVec4(0.3f, 0.3f, 0.3f, 1.0f));

        // ��ʾ������ͼ
        sceneView->display();
        inspectorView->display();
        projectView->display();
	}

    /**
     * ��ʾ������
     * TODO: ʵ�ֹ���������
     */
    void Manager::displayToolBar() {
        // TODO: ʵ�ֹ�����
    }
}