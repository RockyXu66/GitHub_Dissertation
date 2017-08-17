// Std. Includes
#include <string>

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// GL includes
#include "utils/Shader.h"
#include "utils/Camera.h"
#include "utils/Model.h"

// GLM Mathemtics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "PCA.hpp"
#include"Eigen/Eigen"


// Other Libs
#include <SOIL.h>
#include <iostream>

#include "pngwriter.h"

using namespace std;
using namespace glm;
using namespace Eigen;

// Properties
GLuint screenWidth = 800, screenHeight = 600;
int nPixels = 3*screenWidth*screenHeight;

// Function prototypes
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void Do_Movement();
GLuint loadTexture(GLchar* path);
GLuint loadNewTexture(unsigned char* newImage, int width, int height);

// Camera
Camera camera(vec3(0.0f, 1.0f, 5.0f));
bool keys[1024];
GLfloat lastX = 400, lastY = 300;
bool firstMouse = true;

GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

// Global variables
GLuint floorTexture, testTexture;
GLuint planeVAO;

//bool save_screenshot(string filename, int w, int h)
//{
//    //This prevents the images getting padded
//    // when the width multiplied by 3 is not a multiple of 4
//    glPixelStorei(GL_PACK_ALIGNMENT, 1);
//    
//    int nSize = w*h*3;
//    // First let's create our buffer, 3 channels per Pixel
//    char* dataBuffer = (char*)malloc(nSize*sizeof(char));
//    
//    if (!dataBuffer) return false;
//    
//    // Let's fetch them from the backbuffer
//    // We request the pixels in GL_BGR format, thanks to Berzeger for the tip
//    glReadPixels((GLint)0, (GLint)0,
//                 (GLint)w, (GLint)h,
//                 GL_BGR, GL_UNSIGNED_BYTE, dataBuffer);
//    
//    //Now the file creation
//    FILE *filePtr = fopen(filename.c_str(), "wb");
//    if (!filePtr) return false;
//    
//    
//    unsigned char TGAheader[12]={0,0,2,0,0,0,0,0,0,0,0,0};
//    unsigned char header[6] = { w%256,w/256,
//        h%256,h/256,
//        24,0};
//    // We write the headers
//    fwrite(TGAheader,	sizeof(unsigned char),	12,	filePtr);
//    fwrite(header,	sizeof(unsigned char),	6,	filePtr);
//    // And finally our image data
//    fwrite(dataBuffer,	sizeof(GLubyte),	nSize,	filePtr);
//    fclose(filePtr);
//    
//    free(dataBuffer);
//    
//    return true;
//}

// The MAIN function, from here we start our application and run our Game loop
int main()
{
    // Init GLFW
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    
    GLFWwindow* window = glfwCreateWindow(screenWidth, screenHeight, "PCA", nullptr, nullptr); // Windowed
    glfwMakeContextCurrent(window);
    
    // Set the required callback functions
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    
    // Options
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    
    // Initialize GLEW to setup the OpenGL Function pointers
    glewExperimental = GL_TRUE;
    glewInit();
    
    // Define the viewport dimensions
    //glViewport(0, 0, screenWidth, screenHeight);
    
    // Setup some OpenGL options
    glEnable(GL_DEPTH_TEST);
    
    // Setup and compile our shaders
    Shader lightingShader("shaders/basic.vs", "shaders/basic.frag");
    
    GLfloat planeVertices[] = {
        // Positions            // Normals           // Texture Coords
        -1.0f, 1.0f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f, -1.0f,
        1.0f, -1.0f, 0.0f,  0.0f,  1.0f,  0.0f,  -1.0f, 0.0f,
        -1.0f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f, 0.0f,
        
        -1.0f, 1.0f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f, -1.0f,
        1.0f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,  -1.0f, 0.0f,
        1.0f, 1.0f, 0.0f,  0.0f,  1.0f,  0.0f,  -1.0f, -1.0f
    };
    
    //Setup plane VAO
    GLuint planeVBO;
    glGenVertexArrays(1, &planeVAO);
    glGenBuffers(1, &planeVBO);
    glBindVertexArray(planeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, planeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(planeVertices), &planeVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)(6 * sizeof(GLfloat)));
    glBindVertexArray(0);
    
    //Light source
    vec3 lightPos(0.0f, 10.0f, 0.0f);
    
    //Load textures
    floorTexture = loadTexture("textures/container2.png");
    testTexture = loadTexture("textures/image_face/oneColor.png");
    
    
    int image_width, image_height;
    int image_num = 36;
    
//    // Load images
//    vector<unsigned char*> images;
//    for(int i=0; i<image_num; i++){
//        string filename = "textures/new_Model_screenShot/"+to_string(i+1)+".png";
//        unsigned char* image = SOIL_load_image(filename.c_str(), &image_width, &image_height, 0, SOIL_LOAD_RGB);
//        images.push_back(image);
//    }
//
//    int testTime1 = glfwGetTime();
//    
//    // PCA
//    PCA pca;
//    unsigned char* newImage = pca.CalculateEigen(images, image_width, image_height);
//
//    testTexture = loadNewTexture(newImage, image_width, image_height);
//    
//    int testTime2 = glfwGetTime();
//    int duration = testTime2 - testTime1;
//    cout<<"PCA duration time: "<<endl;
//    cout<<duration<<" s"<<endl;
//    cout<<float(duration)/60.0f<<" min"<<endl;
    
/*
//    unsigned char* pixels = SOIL_load_image("textures/image_face/eye_1.png", &image_width, &image_height, 0, SOIL_LOAD_RGBA);
//    vector<unsigned char*> images;
//    images.push_back(pixels);
    
//    vector<int> r_values;
//    vector<int> g_values;
//    vector<int> b_values;
//    vector<int> alpha_values;
    
//    for(int i=0; i<image_width*image_height; i++){
//        int r = images[0][i*4];
//        int g = image[i*4+1];
//        int b = image[i*4+2];
//        int alpha = image[i*4+3];
//        
//        r_values.push_back(r);
//        g_values.push_back(g);
//        b_values.push_back(b);
//        alpha_values.push_back(alpha);
//    }
    
//    pca.CalculateEigen(unsigned char* image);
    
//    for(int i=0; i<r_values.size(); i++){
//        cout<<r_values[i]<<" ";
//    }
//    cout<<endl;
//    cout<<r_values.size();
    
    //    unsigned char (*pixels)[256][256][3] = (unsigned char (*)[256][256][3]) pixels_;
//    int testValue = pixels[0][10];
//    
//    int y,x,c;
//    for (y = 0; y < 20; y++){
//        for (x = 0; x < 20; x++) {
//            int r = pixels[0][x][y];
//            int g = pixels[1][x][y];
//            int b = pixels[2][x][y];
//            int a = pixels[3][x][y];
//            
//            int cc = 0;
//            // etc
//        }
//    }

    // 155 108 78 255
//       108 78 155 0
//       78 155 108 0
//       155 108 78 0
//       108 78 155 0
//       78 155 108 0
    
//    for (int i = 0; i < image_width; i++)
//    {
//        for (int j = 0; j < image_height; j++)
//        {
////            if (pixels_[i + j * image_width] == 1)
////            {
//                //Do Something With Data
//                unsigned char r = pixels_[(i + j * image_width) * 3 + 0];
//                unsigned char g = pixels_[(i + j * image_width) * 3 + 1];
//                unsigned char b = pixels_[(i + j * image_width) * 3 + 2];
//                
//                //...
////            }
//        }
//    }
*/
    
    int count_record = 0;
    
    // Game loop
    while(!glfwWindowShouldClose(window))
    {
//        save_screenshot(record.png, screenWidth, screenHeight);
        GLfloat* pixels = new GLfloat[nPixels*4];
        glReadPixels(0.0, 0.0, screenWidth*2, screenHeight*2,GL_RGB, GL_FLOAT, pixels);
        pngwriter PNG(screenWidth*2, screenHeight*2, 1.0, "record.png");
        size_t x = 1;
        size_t y = 1;
        double R, G, B;
        for(size_t i=0; i<nPixels*4; i++) // "i" is the index for array "pixels"
        {
            switch(i%3)
            {
                case 2:
                    B = static_cast<double>(pixels[i]); break;
                case 1:
                    G = static_cast<double>(pixels[i]); break;
                case 0:
                    R = static_cast<double>(pixels[i]);
                    PNG.plot(x, y, R, G, B);     // set pixel to position (x, y)
                    if( x == screenWidth*2 )             // Move to the next row of image
                    {
                        x=1;
                        y++;
                    }
                    else                       // To the next pixel
                    { x++; }
                    break;
            }
        }
        PNG.close();
        
        // Set frame time
        GLfloat currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        
        GLfloat time = glfwGetTime();
        
        // Check and call events
        glfwPollEvents();
        Do_Movement();
        
        // Clear the colorbuffer
        glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        
        // Transformation matrices
        mat4 projection = perspective(camera.Zoom, (float)screenWidth/(float)screenHeight, 0.1f, 100.0f);
        mat4 view = camera.GetViewMatrix();
        mat4 model;
        model = rotate(model, radians(180.0f), vec3(0.0f, 1.0f, 0.0));
        
        //Draw objects
        lightingShader.Use();
        glUniformMatrix4fv(glGetUniformLocation(lightingShader.Program, "model"), 1, GL_FALSE, value_ptr(model));
        glUniformMatrix4fv(glGetUniformLocation(lightingShader.Program, "view"), 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(lightingShader.Program, "projection"),1,GL_FALSE,value_ptr(projection));
        glBindVertexArray(planeVAO);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, floorTexture);
        glUniform1i(glGetUniformLocation(lightingShader.Program, "floorTexture"), 0);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, testTexture);
        glUniform1i(glGetUniformLocation(lightingShader.Program, "testTexture"), 1);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
        
        // Swap the buffers
        glfwSwapBuffers(window);
    }
    
    glfwTerminate();
    return 0;
}

// This function loads a texture from file. Note: texture loading functions like these are usually
// managed by a 'Resource Manager' that manages all resources (like textures, models, audio).
// For learning purposes we'll just define it as a utility function.
GLuint loadTexture(GLchar* path)
{
    // Generate texture ID and load texture data
    GLuint textureID;
    glGenTextures(1, &textureID);
    int width,height;
    unsigned char* image = SOIL_load_image(path, &width, &height, 0, SOIL_LOAD_RGB);
    
    // Assign texture to ID
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
    glGenerateMipmap(GL_TEXTURE_2D);
    
    // Parameters
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
    SOIL_free_image_data(image);
    return textureID;
}

GLuint loadNewTexture(unsigned char* newImage, int width, int height){
    GLuint textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, newImage);
    glGenerateMipmap(GL_TEXTURE_2D);
    
    // Parameters
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
//    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
    SOIL_free_image_data(newImage);
    return textureID;
}

#pragma region "User input"

// Moves/alters the camera positions based on user input
void Do_Movement()
{
    // Camera controls
    if(keys[GLFW_KEY_W])
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if(keys[GLFW_KEY_S])
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if(keys[GLFW_KEY_A])
        camera.ProcessKeyboard(LEFT, deltaTime);
    if(keys[GLFW_KEY_D])
        camera.ProcessKeyboard(RIGHT, deltaTime);
}

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
    
    if(action == GLFW_PRESS)
        keys[key] = true;
    else if(action == GLFW_RELEASE)
        keys[key] = false;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if(firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }
    
    GLfloat xoffset = xpos - lastX;
    GLfloat yoffset = lastY - ypos;
    
    lastX = xpos;
    lastY = ypos;
    
    camera.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}

#pragma endregion
