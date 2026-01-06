#pragma once
#ifndef __DEBUG_CAPTURE_H__
#define __DEBUG_CAPTURE_H__

#include <string>
#include <fstream>
#include <vector>
#include <glad/glad.h>
#include <glm/glm.hpp>

namespace Glb {

class DebugCapture {
public:
    static DebugCapture& getInstance() {
        static DebugCapture instance;
        return instance;
    }

    void startSession();           // 创建新的 run_YYYY-MM-DD_HHMM 目录
    void endSession();             // 结束会话，生成 clip.mp4

    void captureFrame(GLuint textureId, int frameNum);  // 捕获帧
    void saveConfig(const std::string& methodName);     // 保存当前参数
    void saveCamera();                                   // 保存相机状态
    void log(const std::string& msg);                   // 写日志

    bool isCapturing() const { return capturing; }
    void setCapturing(bool v) { capturing = v; }

    std::string getSessionPath() const { return sessionPath; }

private:
    DebugCapture() : capturing(false), frameCount(0) {}
    DebugCapture(const DebugCapture&) = delete;
    DebugCapture& operator=(const DebugCapture&) = delete;

    bool capturing;
    int frameCount;
    std::string sessionPath;
    std::ofstream logFile;

    void savePNG(const std::string& path, const std::vector<unsigned char>& pixels, int w, int h);
    std::string getCurrentTimestamp();
};

}

#endif
