#pragma once
#ifndef __SCENE_VIEW_H__
#define __SCENE_VIEW_H__

#include <iostream>

#include "glad/glad.h"
#include "glfw3.h"

#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Configure.h"
#include "Component.h"

#include "Global.h"

#include "Manager.h"
#include "DebugCapture.h"
#include <thread>


namespace FluidSimulation {
	class SceneView {

	public:
		GLFWwindow* window;
		ImVec2 pos;

		// std::thread renderThread;

		GLuint texture;

		// �������Ƿ�ͣ������Ⱦ����
		bool isMouseHovering = false;
		// ����϶�״̬
		bool isLeftMouseDragging = false;
		bool isRightMouseDragging = false;
		double lastMouseX = 0.0;
		double lastMouseY = 0.0;
		double mouseX, mouseY;

		int captureFrameNum = 0;

	public:
		SceneView();
		SceneView(GLFWwindow* window);
		void display();
	};
}

#endif