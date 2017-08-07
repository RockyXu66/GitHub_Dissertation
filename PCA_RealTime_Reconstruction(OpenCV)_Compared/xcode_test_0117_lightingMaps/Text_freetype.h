//
//  Text_freetype.h
//  armIK_0221
//
//  Created by Yinghan Xu on 10/03/2017.
//  Copyright © 2017 Yinghan Xu. All rights reserved.
//

#ifndef Text_freetype_h
#define Text_freetype_h

// Std. Includes
#include <string>

// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

// GLFW
#include <GLFW/glfw3.h>

// GL includes
#include "Shader.h"
#include "Camera.h"
#include "Model.h"

// GLM Mathemtics
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

// Other Libs
#include <SOIL.h>
#include "Skeleton.h"
#include "Bone.h"
#include "IKAnalytical.h"
#include "IKCCD.h"
#include "Interpolation.h"
#include "AntTweakBar.h"

// FreeType
#include <ft2build.h>
#include FT_FREETYPE_H


#include <iostream>

using namespace glm;

void Text_freetype(){

// Compile and setup the shader
Shader shader("shaders/text.vs", "shaders/text.frag");
glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(screenWidth), 0.0f, static_cast<GLfloat>(screenHeight));
shader.Use();
glUniformMatrix4fv(glGetUniformLocation(shader.Program, "projection"), 1, GL_FALSE, glm::value_ptr(projection));

// FreeType
FT_Library ft;
// All functions return a value different than 0 whenever an error occurred
if (FT_Init_FreeType(&ft))
std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;

// Load font as face
FT_Face face;
if (FT_New_Face(ft, "fonts/arial.ttf", 0, &face))
std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;

// Set size to load glyphs as
FT_Set_Pixel_Sizes(face, 0, 48);

// Disable byte-alignment restriction
glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

// Load first 128 characters of ASCII set
for (GLubyte c = 0; c < 128; c++)
{
    // Load character glyph
    if (FT_Load_Char(face, c, FT_LOAD_RENDER))
    {
        std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
        continue;
    }
    // Generate texture
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(
                 GL_TEXTURE_2D,
                 0,
                 GL_RED,
                 face->glyph->bitmap.width,
                 face->glyph->bitmap.rows,
                 0,
                 GL_RED,
                 GL_UNSIGNED_BYTE,
                 face->glyph->bitmap.buffer
                 );
    // Set texture options
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // Now store character for later use
    Character character = {
        texture,
        glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
        glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
        face->glyph->advance.x
    };
    Characters.insert(std::pair<GLchar, Character>(c, character));
}
glBindTexture(GL_TEXTURE_2D, 0);
// Destroy FreeType once we're finished
FT_Done_Face(face);
FT_Done_FreeType(ft);

// Configure VAO/VBO for texture quads
glGenVertexArrays(1, &VAOText);
glGenBuffers(1, &VBOText);
glBindVertexArray(VAOText);
glBindBuffer(GL_ARRAY_BUFFER, VBOText);
glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
glEnableVertexAttribArray(0);
glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);
glBindBuffer(GL_ARRAY_BUFFER, 0);
glBindVertexArray(0);
}
#endif /* Text_freetype_h */
