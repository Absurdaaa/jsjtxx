#pragma once
#ifndef __MANAGER_H__
#define __MANAGER_H__

// OpenGL���ͷ�ļ�
#include "glad/glad.h"
#include "glfw3.h"

// ImGui���ͷ�ļ�
#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// UI��ͼͷ�ļ�
#include "SceneView.h"
#include "InspectorView.h"
#include "ProjectView.h"

// ����ģ�����ͷ�ļ�
#include "Lagrangian2dComponent.h"
#include "Lagrangian2dFountainComponent.h"
#include "Eulerian2dComponent.h"
#include "PICComponent.h"
#include "Lagrangian3dComponent.h"
#include "Eulerian3dComponent.h"
#include "PIC3dComponent.h"

#include <vector>

namespace FluidSimulation
{
	class SceneView;
	class InspectorView;
	class ProjectView;

	// UI��������(����ģʽ)
	// �����������UI��ͼ������ģ�����
	class Manager {
	public:
		// ��ȡ����ʵ��
		static Manager& getInstance() {
			static Manager instance;
			return instance;
		}

		void init(GLFWwindow* window);      // ��ʼ��������
		void displayViews();                // ��ʾ������ͼ
		void displayToolBar();              // ��ʾ������

		// Getter/Setter����
		SceneView* getSceneView() const { return sceneView; };
		InspectorView* getInspectorView() const { return inspectorView; };
		ProjectView* getProjectView() const { return projectView; };
		GLFWwindow* getWindow() const { return window; };
		Glb::Component* getMethod() const { return currentMethod; };
		void setMethod(Glb::Component* method) { currentMethod = method; };

	private:
		// ˽�й��캯��(����ģʽ)
		Manager() {
			window = NULL;
			sceneView = NULL;
			inspectorView = NULL;
			projectView = NULL;
			currentMethod = NULL;
		};

		// ��ֹ�����͸�ֵ(����ģʽ)
		Manager(const Manager&) = delete;
		Manager& operator=(const Manager&) = delete;
		
		GLFWwindow* window;                 // GLFW����
		SceneView* sceneView;               // ������ͼ
		InspectorView* inspectorView;       // ��������ͼ
		ProjectView* projectView;           // ��Ŀ��ͼ

		Glb::Component* currentMethod;      // ��ǰѡ���ģ�ⷽ��
	};
}

#endif