// Std. Includes
#include <string>

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

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

using namespace std;
using namespace glm;
using namespace Eigen;
using namespace cv;

// Properties
GLuint screenWidth = 900, screenHeight = 750;

// Function prototypes
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void Do_Movement();
GLuint loadTexture(GLchar* path);
GLuint loadNewTexture(unsigned char* newImage, int width, int height);

// Camera
Camera camera(vec3(0.0f, 0.0f, 2.0f));
bool keys[1024];
GLfloat lastX = 900, lastY = 750;
bool firstMouse = true;

GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

//Global variables
GLuint floorTexture, testTexture;
GLuint planeVAO;
GLuint imageIndex = 10;

// The MAIN function, from here we start our application and run our Game loop
int main(int argc, const char *argv[])
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

    char* path = new char[oneImagePath.length() + 1];
    strcpy(path, oneImagePath.c_str());
    // do stuff
    delete [] path;
    
    //Load textures
    floorTexture = loadTexture(path);
    testTexture = loadTexture(path);
    
    // PCA
    PCA_ pca;
    GLfloat testTime1 = glfwGetTime();
    pca.load();
    GLfloat testTime2 = glfwGetTime();
    GLfloat duration = testTime2 - testTime1;
    int image_width, image_height;
    unsigned char* newImage = SOIL_load_image(oneImagePath.c_str(), &image_width, &image_height, 0, SOIL_LOAD_RGB);
    
    Mat bgr[3];
    split(imread(oneImagePath), bgr);
    
    int cell_num = (image_width*image_height)/(cell_dimension*cell_dimension);
    int x = image_height/cell_dimension;
    int y = image_width/cell_dimension;
    
    GLfloat totalTime = 0;
    int count = 0;
    // Game loop
    while(!glfwWindowShouldClose(window))
    {
        // Set frame time
        GLfloat currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        
        if(deltaTime<0.25){
            count++;
            totalTime += deltaTime;
        }
        
//        cout<<"DeltaTime: "<<deltaTime<<" FPS: "<<1.0f/deltaTime<<endl;
//        cout<<imageIndex<<endl;
        
        // Check and call events
        glfwPollEvents();
        Do_Movement();
        
        // Clear the colorbuffer
        glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        newImage = pca.reconstruct(imageIndex, newImage, bgr, cell_num, x, y);
        testTexture = loadNewTexture(newImage, image_width, image_height);
        
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
    cout<<"====================="<<endl;
    cout<<"DeltaTime: "<<totalTime/float(count)<<" FPS: "<<1.0f/(totalTime/float(count))<<endl;
    string name = "head"+to_string(image_width)+"_cells"+to_string(cell_dimension)+"_"+to_string(num_components)+".txt";
    if(FPS_record){
        FileStorage nfs("/Users/yinghanxu/Study/Dissertation_ResultData/SSIM&FPS/"+name, FileStorage::WRITE);
        nfs << "FPS" << float(1.0f/(totalTime/float(count)));
        nfs.release();
        cout<<endl<<"====Save average FPS to the file===="<<endl<<endl;
    }
    
//    float test;
//    FileStorage fs("/Users/yinghanxu/Study/Dissertation_ResultData/SSIM&FPS/"+name,FileStorage::READ);
//    fs["FPS"] >> test;
//    fs.release();
//    cout<<test<<endl;
    
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
//    SOIL_free_image_data(newImage);   // Keep this
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
    if(keys[GLFW_KEY_RIGHT]){
        imageIndex --;
        if(imageIndex==800){
            imageIndex = 700;
        }else if(imageIndex%100==0){
            imageIndex += 100;
        }
    }
    if(keys[GLFW_KEY_LEFT]){
        imageIndex ++;
        if(imageIndex==900){
            imageIndex = 800;
        }else if((imageIndex-1)%100==0){
            imageIndex -= 100;
        }
    }
    if(keys[GLFW_KEY_UP]){
        if(imageIndex<400||(imageIndex>500&&imageIndex<800)){
            imageIndex += 100;
        }
        if(imageIndex>800){
            imageIndex -= 800;
        }
    }
    if(keys[GLFW_KEY_DOWN]){
        if((imageIndex<900&&imageIndex>=600)||(imageIndex<=500&&imageIndex>=100)){
            imageIndex -= 100;
        }
        if(imageIndex<100){
            imageIndex += 800;
        }
    }
}

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
    
    if(action == GLFW_PRESS){
//        switch (key) {
//                
//            case GLFW_KEY_LEFT:
//                if(imageIndex==0){
//                    cout<<"imageIndex: "<<imageIndex;
//                    imageIndex = 119;
//                }else{
//                    imageIndex--;
//                }
//                break;
//            case GLFW_KEY_RIGHT:
//                imageIndex ++;
//                if(imageIndex > 119){
//                    imageIndex = 0;
//                }
//                break;
//                
//            default:
//                break;
//        }
        keys[key] = true;
    }
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
