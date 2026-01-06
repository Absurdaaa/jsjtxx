#include "DebugCapture.h"
#include "Configure.h"
#include "Camera.h"
#include "Global.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <direct.h>  // Windows: _mkdir

namespace Glb {

// 递归创建目录 (Windows)
static void createDirectories(const std::string& path) {
    size_t pos = 0;
    while ((pos = path.find_first_of("/\\", pos + 1)) != std::string::npos) {
        _mkdir(path.substr(0, pos).c_str());
    }
    _mkdir(path.c_str());
}

std::string DebugCapture::getCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm tm;
    localtime_s(&tm, &time);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H%M");
    return oss.str();
}

void DebugCapture::startSession() {
    sessionPath = "run_" + getCurrentTimestamp();
    createDirectories(sessionPath + "/frames");

    logFile.open(sessionPath + "/log.txt");
    logFile << "Session started: " << getCurrentTimestamp() << "\n";
    logFile << "Resolution: " << imageWidth << "x" << imageHeight << "\n\n";

    frameCount = 0;
    capturing = true;
    log("Capture session started");
}

void DebugCapture::endSession() {
    if (!capturing) return;

    log("Session ended, total frames: " + std::to_string(frameCount));
    logFile.close();
    capturing = false;

    // 生成视频 (需要 ffmpeg 在 PATH 中)
    std::string cmd = "ffmpeg -y -framerate 30 -i \"" + sessionPath + "/frames/%06d_final.png\" "
                      "-c:v libx264 -pix_fmt yuv420p \"" + sessionPath + "/clip.mp4\" 2>nul";
    system(cmd.c_str());
}

void DebugCapture::captureFrame(GLuint textureId, int frameNum) {
    if (!capturing || textureId == 0) return;

    std::vector<unsigned char> pixels(imageWidth * imageHeight * 4);

    glBindTexture(GL_TEXTURE_2D, textureId);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    glBindTexture(GL_TEXTURE_2D, 0);

    // 翻转Y轴 (OpenGL坐标系)
    std::vector<unsigned char> flipped(pixels.size());
    for (int y = 0; y < imageHeight; y++) {
        memcpy(&flipped[y * imageWidth * 4],
               &pixels[(imageHeight - 1 - y) * imageWidth * 4],
               imageWidth * 4);
    }

    std::ostringstream filename;
    filename << sessionPath << "/frames/" << std::setfill('0') << std::setw(6) << frameNum << "_final.png";
    savePNG(filename.str(), flipped, imageWidth, imageHeight);

    frameCount++;
}

void DebugCapture::savePNG(const std::string& path, const std::vector<unsigned char>& pixels, int w, int h) {
    stbi_write_png(path.c_str(), w, h, 4, pixels.data(), w * 4);
}

void DebugCapture::saveConfig(const std::string& methodName) {
    if (sessionPath.empty()) return;

    std::ofstream f(sessionPath + "/config.json");
    f << "{\n";
    f << "  \"method\": \"" << methodName << "\",\n";
    f << "  \"resolution\": [" << imageWidth << ", " << imageHeight << "],\n";

    if (methodName.find("Eulerian2d") != std::string::npos) {
        f << "  \"dt\": " << Eulerian2dPara::dt << ",\n";
        f << "  \"contrast\": " << Eulerian2dPara::contrast << ",\n";
        f << "  \"gridNum\": " << Eulerian2dPara::gridNum << ",\n";
        f << "  \"boussinesqAlpha\": " << Eulerian2dPara::boussinesqAlpha << ",\n";
        f << "  \"boussinesqBeta\": " << Eulerian2dPara::boussinesqBeta << "\n";
    }
    else if (methodName.find("Eulerian3d") != std::string::npos) {
        f << "  \"dt\": " << Eulerian3dPara::dt << ",\n";
        f << "  \"contrast\": " << Eulerian3dPara::contrast << ",\n";
        f << "  \"gridNum\": [" << Eulerian3dPara::gridNumX << "," << Eulerian3dPara::gridNumY << "," << Eulerian3dPara::gridNumZ << "],\n";
        f << "  \"boussinesqAlpha\": " << Eulerian3dPara::boussinesqAlpha << ",\n";
        f << "  \"boussinesqBeta\": " << Eulerian3dPara::boussinesqBeta << "\n";
    }
    else if (methodName.find("Lagrangian2d") != std::string::npos) {
        f << "  \"dt\": " << Lagrangian2dPara::dt << ",\n";
        f << "  \"substep\": " << Lagrangian2dPara::substep << ",\n";
        f << "  \"viscosity\": " << Lagrangian2dPara::viscosity << ",\n";
        f << "  \"stiffness\": " << Lagrangian2dPara::stiffness << ",\n";
        f << "  \"density\": " << Lagrangian2dPara::density << "\n";
    }
    else if (methodName.find("Lagrangian3d") != std::string::npos) {
        f << "  \"dt\": " << Lagrangian3dPara::dt << ",\n";
        f << "  \"substep\": " << Lagrangian3dPara::substep << ",\n";
        f << "  \"viscosity\": " << Lagrangian3dPara::viscosity << ",\n";
        f << "  \"stiffness\": " << Lagrangian3dPara::stiffness << ",\n";
        f << "  \"xsph_c\": " << Lagrangian3dPara::xsph_c << "\n";
    }
    else if (methodName.find("PIC") != std::string::npos && methodName.find("2d") != std::string::npos) {
        f << "  \"dt\": " << PIC2dPara::dt << ",\n";
        f << "  \"contrast\": " << PIC2dPara::contrast << ",\n";
        f << "  \"particlesPerStep\": " << PIC2dPara::particlesPerStep << "\n";
    }
    else if (methodName.find("PIC") != std::string::npos && methodName.find("3d") != std::string::npos) {
        f << "  \"dt\": " << PIC3dPara::dt << ",\n";
        f << "  \"contrast\": " << PIC3dPara::contrast << ",\n";
        f << "  \"particlesPerStep\": " << PIC3dPara::particlesPerStep << "\n";
    }
    else {
        f << "  \"note\": \"unknown method\"\n";
    }

    f << "}\n";
    f.close();
    log("Config saved");
}

void DebugCapture::saveCamera() {
    if (sessionPath.empty()) return;

    auto& cam = Camera::getInstance();
    std::ofstream f(sessionPath + "/camera.json");
    f << "{\n";
    f << "  \"position\": [" << cam.mPosition.x << ", " << cam.mPosition.y << ", " << cam.mPosition.z << "],\n";
    f << "  \"yaw\": " << cam.mYaw << ",\n";
    f << "  \"pitch\": " << cam.mPitch << ",\n";
    f << "  \"fov\": " << cam.fovyDeg << ",\n";
    f << "  \"near\": " << cam.nearPlane << ",\n";
    f << "  \"far\": " << cam.farPlane << "\n";
    f << "}\n";
    f.close();
    log("Camera saved");
}

void DebugCapture::log(const std::string& msg) {
    if (!logFile.is_open()) return;

    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm tm;
    localtime_s(&tm, &time);

    logFile << "[" << std::put_time(&tm, "%H:%M:%S") << "] " << msg << "\n";
    logFile.flush();
}

}
