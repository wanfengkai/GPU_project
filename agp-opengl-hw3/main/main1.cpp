#include "common.hpp"
#include "util.hpp"

using namespace std;
using namespace glm;
using namespace agp;

GLuint g_default_vao = 0;
GLuint shaderProgram;
GLuint MLoc, VPLoc, MVPLoc, inputColorLoc;

const unsigned int width = 1280;
const unsigned int height = 720;
const int NUM_PARTICLES = 100;

float fov = 45.0f;
float angle1 = 0.0f;
int dt = 1;

struct Particle
{
    glm::vec3 p;
    glm::vec3 v;
};

Particle Pa[NUM_PARTICLES];


void initial(Particle *Pa){
   
    for (int i = 0; i<NUM_PARTICLES; i++)  
    {
        
        Pa[i].p.x =   5 * i / NUM_PARTICLES + 5 * rand() / float(RAND_MAX) ;
        Pa[i].p.y =   i / NUM_PARTICLES + 5 * rand() / float(RAND_MAX) ;
        Pa[i].p.z =   i / NUM_PARTICLES + 5 * rand() / float(RAND_MAX) ;
    }
}
/*
void uPaate(Particle *Pa){
   
    for (int i = 0; i<NUM_PARTICLES; i++)  
    {
        if(i % 4 == 0){
        Pa[i].v.x =  0.0001f* rand() / float(RAND_MAX) ;
        Pa[i].v.y =  0.0001f* rand() / float(RAND_MAX) ;
        Pa[i].v.z =  0.0f ;
    }
        if(i % 4 == 1){
        Pa[i].v.x =  -0.0001f* rand() / float(RAND_MAX) ;
        Pa[i].v.y =  0.0001f* rand() / float(RAND_MAX) ;
        Pa[i].v.z =  0.0f ;
    }
        if(i % 4 == 2){
        Pa[i].v.x =  0.0001f* rand() / float(RAND_MAX) ;
        Pa[i].v.y =  -0.0001f* rand() / float(RAND_MAX) ;
        Pa[i].v.z =  0.0f ;
    }
        if(i % 4 == 3){
        Pa[i].v.x =  -0.0001f* rand() / float(RAND_MAX) ;
        Pa[i].v.y =  -0.0001f* rand() / float(RAND_MAX) ;
        Pa[i].v.z =  0.0f ;
    }
        Pa[i].p.x += Pa[i].v.x * dt;
        Pa[i].p.y += Pa[i].v.y * dt;
        Pa[i].p.z += Pa[i].v.z * dt; 
    }
}
    */
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
}

void release()
{
    // Release the default VAO
    glDeleteVertexArrays(1, &g_default_vao);
    
    // Do not forget to release any memory allocation here!
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------

void display()
{
    // Clear the screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //printf("FreeGLUT triggered the display() callback!\n");
    
    // Your rendering code must be here! Do not forget to swap the
    // front and back buffers, and to force a redisplay to keep the
    // render loop running. This functionality is available within
    // FreeGLUT as well, check the assignment for more information.
    
    // Important note: The following function flushes the rendering
    // queue, but this is only for single-buffered rendering. You
    // must replace this function following the previous indications.
    //glFlush();
    


    mat4 M,V,P,MVP,VP; 
    mat4 rotationMat(1);  
    rotationMat = rotate(rotationMat, angle1, vec3(0.0, 1.0, 0.0));
    vec3 vec = vec3(rotationMat * vec4(0.0f, 0.0f, 9.0f, 1.0f)); 
    
    V = lookAt(vec, vec3(0.0f, 0.0f, 0.0f),  vec3(0.0f, 1.0f, 0.0f));
    P = perspective(radians(fov), (float)width/(float)height, 0.1f, 100.0f);
    VP = P * V;
    
    glUniformMatrix4fv(VPLoc, 1, GL_FALSE, glm::value_ptr(VP)); 
    


    for (int i = 0; i<NUM_PARTICLES; i++)  
    {  
        vec4 inputColor = vec4(0.5f - ((float)i/NUM_PARTICLES)/2.f, ((float)i/NUM_PARTICLES), 0.5f + ((float)i/NUM_PARTICLES)/2.0f, 0.5f);
        
        //srand(i);  
        //float r= 10 * rand() / float(RAND_MAX);
        
        //update(Pa);
        launchkernel(Particle *Pa)

        M = mat4();// By default, this is identity matrix
        M = translate(M, vec3(Pa[i].p.x, Pa[i].p.y, Pa[i].p.z));
        //MVP = P*V*M;

        // pass them to the shaders

        glUniformMatrix4fv(MLoc, 1, GL_FALSE, glm::value_ptr(M));
        glUniform4fv(inputColorLoc, 1, glm::value_ptr(inputColor));

        glutSolidSphere(r,20,20);
        glutWireSphere(r,20,20);
    }
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
         angle1 += 0.1f;
    }  
 
    if(key==GLUT_KEY_RIGHT)  
    {  
         angle1 -= 0.1f;
    }  
}  




int main(int argc, char **argv)
{

    initial(Pa);
    
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

