#version 330 core
layout (location = 0) in vec3 aPos;
 
uniform mat4 M;
uniform mat4 VP;

out vec4 gl_Position;

void main(){
	
	mat4 MVP = VP * M;
	gl_Position = MVP * vec4(aPos, 1);

}

