#version 330 core

uniform sampler2D texture_diffuse;

in vec2 TexCoords;

out vec4 color;

void main()
{           
    color = vec4(texture(texture_diffuse, TexCoords));
}



