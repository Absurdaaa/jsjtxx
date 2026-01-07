#pragma once
#ifndef __INSPECTOR_VIEW_H__
#define __INSPECTOR_VIEW_H__

// OpenGL���ͷ�ļ�
#include "glad/glad.h"
#include "glfw3.h"

// ImGui���ͷ�ļ�
#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// ��Ŀ���ͷ�ļ�
#include "Configure.h"
#include "Manager.h"
#include "Logger.h"
#include "DebugCapture.h"

// InspectorView.cpp 直接访问 Glb::Camera 的成员，确保此处可见完整定义。
#include "Camera.h"

#include <iostream>
#include <string>

namespace FluidSimulation {
	// ��������ͼ��
	// ������ʾ�ͱ༭��ǰ����ģ�ⷽ���Ĳ���
	class InspectorView {
	private:
		GLFWwindow* window;	// GLFW����
		ImVec2 pos;			// ��ͼλ��

	public:
		int showID;			// �Ƿ���ʾ���ID

		InspectorView();
		InspectorView(GLFWwindow* window);
		void display();		// ��ʾ�����༭���
	};
}

#endif