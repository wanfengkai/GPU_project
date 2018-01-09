// #include "common.hpp"
// #include "util.hpp"
// #include "cuda_kernel.cu"

#include "part.h"
#include <time.h>
#include <stdio.h>
#include <string.h>

#define D 376.78

// # include <cstdlib>

// wait for teacher

//-------------------------------------------------------------------//


using namespace std;
using namespace glm;
using namespace agp;

GLuint g_default_vao = 0;
GLuint shaderProgram;
GLuint MLoc, VPLoc, MVPLoc, inputColorLoc;

const unsigned int width = 1280;
const unsigned int height = 720;


float fov = 45.0f;
float rotangle = 0.0f;
int dt = 1;
vec4 red = vec4(0.5f, 0.0f, 0.0f, 0.5f);      // color red, silicate particle
vec4 yellow = vec4(0.5f, 0.5f, 0.0f, 0.5f);      // color yellow, iron particle
vec4 inputColor = vec4(0.0f, 0.0f, 0.0f, 0.0f);

bool init_display = false;
struct Particle *cpuP = (struct Particle *)malloc(P_NUM*sizeof(struct Particle));
// struct vecfloat3 *cpuA = (struct vecfloat3 *)malloc(P_NUM*sizeof(struct vecfloat3));


void init()
{
    // Generate and bind the default VAO
    glGenVertexArrays(1, &g_default_vao);
    glBindVertexArray(g_default_vao);
   
    // Set the background color (RGBA)
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    
    // Your OpenGL settings, such as alpha, depth and others, should be
    // defined here! For the assignment, we only ask you to enable the
    // alpha channel.
    shaderProgram = util::loadShaders("vertex.glsl","fragment.glsl");
    glUseProgram(shaderProgram);

    // retrieve the matrix uniform locations
    MLoc = glGetUniformLocation(shaderProgram, "M"); 
    VPLoc = glGetUniformLocation(shaderProgram, "VP");
    inputColorLoc = glGetUniformLocation(shaderProgram, "inputColor");
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	init_display = true;
}


void draw_particle()
{
    mat4 M,V,P,MVP,VP; 
    //matrix V
    mat4 rotationMat(1);  
    rotationMat = rotate(rotationMat, rotangle, vec3(0.0, 1.0, 0.0));
    vec3 vec = vec3(rotationMat * vec4(0.0f, 0.0f, 50.0f, 1.0f));   
    V = lookAt(vec, vec3(0.0f, 0.0f, 0.0f),  vec3(0.0f, 1.0f, 0.0f));
    //matrix P
    P = perspective(radians(fov), (float)width/(float)height, 0.1f, 10000.0f);
    //matrix VP
    VP = P * V;  
    glUniformMatrix4fv(VPLoc, 1, GL_FALSE, glm::value_ptr(VP)); 

    for (int i = 0; i<P_NUM; i++)  
    {  

		float r = D/20000.0;

        // By default, this is identity matrix
        M = mat4();

        M = translate(M, vec3(cpuP[i].position.x/10000.0, cpuP[i].position.y/10000.0, cpuP[i].position.z/10000.0));
 		
 		if (cpuP[i].p_type == true)
 			inputColor = red;
 		else
 			inputColor = yellow;

        glUniformMatrix4fv(MLoc, 1, GL_FALSE, glm::value_ptr(M));
        glUniform4fv(inputColorLoc, 1, glm::value_ptr(inputColor));

        glutSolidSphere(r,20,20);

    }
//    printf("drawing particles\n");

}

void release()
{
    // Release the default VAO
    glDeleteVertexArrays(1, &g_default_vao);
    
    // Do not forget to release any memory allocation here!
//	cuda_release();
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------

void display()
{
    // Clear the screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

 //   glLoadIdentity();
    
 //   glBindVertexArray(g_default_vao);

    
    //printf("FreeGLUT triggered the display() callback!\n");
    
    // Your rendering code must be here! Do not forget to swap the
    // front and back buffers, and to force a redisplay to keep the
    // render loop running. This functionality is available within
    // FreeGLUT as well, check the assignment for more information.
    
    // Important note: The following function flushes the rendering
    // queue, but this is only for single-buffered rendering. You
    // must replace this function following the previous indications.
    //glFlush();

    if(init_display == true)
    {
        time_t t;
    	time(&t);
		particle_init(t, cpuP);
		printf("particle initialization\n");
        init_display = false;
    }
//    printf("before particle position is: %f, %f, %f\n", cpuP[10].position.x, cpuP[10].position.y, cpuP[10].position.z);
    draw_particle();
    particle_update(cpuP);
//    printf("after particle position is: %f, %f, %f\n", cpuP[10].position.x, cpuP[10].position.y, cpuP[10].position.z);
//   printf("particle update\n");
    // render loop

    glutSwapBuffers();
    glutPostRedisplay();
}


void ProcessNormalKeys(unsigned char key,int x,int y)
{
    if(key==27)//esc
    {
        release();
        exit(0);
    }
    else if(key==43)//+
    {
        fov -= 2.0f;
    }
            
    else if(key==45)//-
    {
        fov += 2.0f;
    }        
    else 
       printf("%c\n",key); 
} 


 
void ProcessSpecialKeys(GLint key,GLint x,GLint y)
{  
      
    if(key==GLUT_KEY_LEFT)  
    {  
         rotangle += 0.1f;
    }  
 
    if(key==GLUT_KEY_RIGHT)  
    {  
         rotangle -= 0.1f;
    }  
}  




int main(int argc, char **argv)
{
    
    // Initialize FreeGLUT and create the window
    glutInit(&argc, argv);
    
    // Setup the window (e.g., size, display mode and so on)
    // glutInitWindowSize( ... );
    // glutInitWindowPosition( ... );
    // glutInitDisplayMode( ... );
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(width, height); 
    
    // Make FreeGLUT return from the main rendering loop
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
            GLUT_ACTION_GLUTMAINLOOP_RETURNS);

    // Create the window and associate the callbacks
    glutCreateWindow("Applied GPU Programming");
    glutDisplayFunc(display);
    // glutIdleFunc( ... );
    // glutReshapeFunc( ... );
    // glutKeyboardFunc( ... );
     glutSpecialFunc(ProcessSpecialKeys);
    // glutMouseFunc( ... );
    // glutMotionFunc( ... );
    glutKeyboardFunc(ProcessNormalKeys);
    

    
    // Init GLAD to be able to access the OpenGL API
    if (!gladLoadGL())
    {
      return GL_INVALID_OPERATION;
    }


    
    // Display OpenGL information
    util::displayOpenGLInfo();
    
    // Initialize the 3D view
    init();
    
    // Launch the main loop for rendering
    glutMainLoop();
    
    // Release all the allocated memory
    release();
    
    return 0;
}

