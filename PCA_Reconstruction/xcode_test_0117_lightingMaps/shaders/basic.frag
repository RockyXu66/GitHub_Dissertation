#version 330 core

uniform sampler2D floorTexture;
uniform sampler2D testTexture;

in vec2 TexCoords;

out vec4 color;

void main()
{           
    vec3 temp = texture(testTexture, TexCoords).rgb;
    color = vec4(temp, 1.0f);
}



