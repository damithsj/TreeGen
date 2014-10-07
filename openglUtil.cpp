#include <stdlib.h> //Needed for "exit" function

//Include OpenGL header files, so that we can use OpenGL
#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "openglUtil.h"

float stepSize = 0.3f;		// Step size
float angleStepSize = 3.0f;		// Angle Step size

float movementX = 0.0f;
float movementY = -50.0f;
float movementZ = -180.0f;

float rotateAngle = -120.0f;
float azimuthAngle = -90.0f; // initial 90 degrees rotation along x azis is needed since the tree growth happens in Z direction
//float rotateAngle = 0.0f;
//float azimuthAngle = 0.0f; // initial 90 degrees rotation along x azis is needed since the tree growth happens in Z direction

//Called when a key is pressed
void handleKeypress(unsigned char key, int x, int y) 
{
	switch (key) {
	case 27: //Escape key
		exit(0); //Exit the program
		break;
	case 'w': //move up
		movementY += stepSize;
		break;
	case 's': //move down
		movementY -= stepSize;
		break;
	case 'a': //move left
		movementX -= stepSize;
		break;
	case 'd': //move right
		movementX += stepSize;
		break;
	case 'r':
		movementZ -= stepSize; //zoom out
		break;
		case 't':
		movementZ += stepSize; //zoom in
		break;

	}

	glutPostRedisplay();
}

void handleSpecialKeyPress(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:		
		azimuthAngle += angleStepSize;
		break;
	case GLUT_KEY_DOWN:		
		azimuthAngle -= angleStepSize;
		break;
	case GLUT_KEY_LEFT:		
		rotateAngle += angleStepSize;
		break;
	case GLUT_KEY_RIGHT:		
		rotateAngle -= angleStepSize;
		break;
	}
	glutPostRedisplay();
}

//Called when the window is resized
void handleResize(int w, int h) 
{
//Tell OpenGL how to convert from coordinates to pixel values
    glViewport(0, 0, w, h);
    
    glMatrixMode(GL_PROJECTION); //Switch to setting the camera perspective
    
    //Set the camera perspective
    glLoadIdentity(); //Reset the camera
    gluPerspective(45.0,                  //The camera angle
                   (double)w / (double)h, //The width-to-height ratio
                   1.0,                   //The near z clipping coordinate
                   200.0);                //The far z clipping coordinate
	
}

//Initializes 3D rendering
void initRendering(int Width, int Height) 
{
   glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
   glColor3f( 0.2, 1.0, 0.2 );					// Greenish color
   glClearDepth(1.0);
   glDepthFunc(GL_LESS);
   glEnable(GL_DEPTH_TEST);
   //glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE );
   glShadeModel(GL_SMOOTH);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
   glMatrixMode(GL_MODELVIEW);   

   glEnable(GL_LIGHTING);
   GLfloat LightAmbient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
   GLfloat LightDiffuse[] = { 0.0f, 1.0f, 0.0f, 1.0f };
   GLfloat LightPosition[] = { 0.0f, 0.0f, 2.0f, 1.0f };
   glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
   glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
   glLightfv(GL_LIGHT1, GL_POSITION,LightPosition);
   glEnable(GL_LIGHT1);
}

void setColors()
{
	glClearColor(1.0,0.5,0.0,0.5);      // set the clear color
	glColor3f(0.0,1.0,0.0);             // set the drawing color
}