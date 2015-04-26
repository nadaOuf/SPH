/** File:		sph_main.cpp
 ** Author:		Dongli Zhang
 ** Contact:	dongli.zhang0129@gmail.com
 **
 ** Copyright (C) Dongli Zhang 2013
 **
 ** This program is free software;  you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 2 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY;  without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 ** the GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program;  if not, write to the Free Software 
 ** Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "sph_header.h"
#include "sph_data.h"
#include "sph_timer.h"
#include "sph_system.h"
#include  <GL\glew.h>
#include  <GL\glut.h>
#include "slef_def.h"


//#pragma comment(lib, "glew32.lib") 

SPHSystem *sph;

Timer *sph_timer;
char *window_title;

GLuint v;
GLuint f;
GLuint p;

void set_shaders()
{
	char *vs=NULL;
	char *fs=NULL;

	vs=(char *)malloc(sizeof(char)*10000);
	fs=(char *)malloc(sizeof(char)*10000);
	memset(vs, 0, sizeof(char)*10000);
	memset(fs, 0, sizeof(char)*10000);

	FILE *fp;
	char c;
	int count;

	fp=fopen("shader/shader.vs", "r");
	count=0;
	while((c=fgetc(fp)) != EOF)
	{
		vs[count]=c;
		count++;
	}
	fclose(fp);

	fp=fopen("shader/shader.fs", "r");
	count=0;
	while((c=fgetc(fp)) != EOF)
	{
		fs[count]=c;
		count++;
	}
	fclose(fp);

	v=glCreateShader(GL_VERTEX_SHADER);
	f=glCreateShader(GL_FRAGMENT_SHADER);

	const char *vv;
	const char *ff;
	vv=vs;
	ff=fs;

	glShaderSource(v, 1, &vv, NULL);
	glShaderSource(f, 1, &ff, NULL);

	int success;

	glCompileShader(v);
	glGetShaderiv(v, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		char info_log[5000];
		glGetShaderInfoLog(v, 5000, NULL, info_log);
		printf("Error in vertex shader compilation!\n");
		printf("Info Log: %s\n", info_log);
	}

	glCompileShader(f);
	glGetShaderiv(f, GL_COMPILE_STATUS, &success);
	if(!success)
	{
		char info_log[5000];
		glGetShaderInfoLog(f, 5000, NULL, info_log);
		printf("Error in fragment shader compilation!\n");
		printf("Info Log: %s\n", info_log);
	}

	p=glCreateProgram();
	glAttachShader(p, v);
	glAttachShader(p, f);
	glLinkProgram(p);
	glUseProgram(p);

	free(vs);
	free(fs);
}

void draw_box(float ox, float oy, float oz, float width, float height, float length)
{
    glLineWidth(1.0f);
	glColor3f(1.0f, 0.0f, 0.0f);

    glBegin(GL_LINES);   
		
        glVertex3f(ox, oy, oz);
        glVertex3f(ox+width, oy, oz);

        glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox, oy, oz);
        glVertex3f(ox, oy, oz+length);

        glVertex3f(ox+width, oy, oz);
        glVertex3f(ox+width, oy+height, oz);

        glVertex3f(ox+width, oy+height, oz);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox, oy, oz+length);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox, oy+height, oz);

        glVertex3f(ox+width, oy, oz);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox, oy, oz+length);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox+width, oy+height, oz);
        glVertex3f(ox+width, oy+height, oz+length);

        glVertex3f(ox+width, oy+height, oz+length);
        glVertex3f(ox+width, oy, oz+length);

        glVertex3f(ox, oy+height, oz+length);
        glVertex3f(ox+width, oy+height, oz+length);

    glEnd();
}

void init_sph_system()
{
	real_world_origin.x=-10.0f;
	real_world_origin.y=-10.0f;
	real_world_origin.z=-10.0f;

	real_world_side.x=20.0f;
	real_world_side.y=20.0f;
	real_world_side.z=20.0f;

	sph=new SPHSystem();
	sph->init_system();

	sph_timer=new Timer();
	window_title=(char *)malloc(sizeof(char)*50);
}

void init()
{
	 glewInit();
	 GLfloat sun_direction[] = { 0.0, 2.0, 1.0, 1.0 };
     GLfloat sun_intensity[] = { 0.7, 0.7, 0.7, 1.0 };
     GLfloat ambient_intensity[] = { 0.3, 0.3, 0.3, 1.0 };

	  glEnable(GL_LIGHTING);              // Set up ambient light.
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_intensity);
     
      glEnable(GL_LIGHT0);                // Set up sunlight.
      glLightfv(GL_LIGHT0, GL_POSITION, sun_direction);
      glLightfv(GL_LIGHT0, GL_DIFFUSE, sun_intensity);
    
	  
      glEnable(GL_COLOR_MATERIAL);        // Configure glColor().
      glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	  GLfloat ambient[]={0.2,0.2,0.2,1.0};
      GLfloat diffuse[]={1.0,0.8,0.0,1.0};
      GLfloat specular[]={1.0,1.0,1.0,1.0};
      GLint shine = 100;
 
    // glMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
    // glMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
    // glMaterialfv(GL_FRONT,GL_SPECULAR,specular);
    // glMateriali(GL_FRONT,GL_SHININESS,shine);
	 ;
	glViewport(0, 0, window_width, window_height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0, (float)window_width/window_height, 10.0f, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -3.0f);
}

void init_ratio()
{
	sim_ratio.x=real_world_side.x/sph->world_size.x;
	sim_ratio.y=real_world_side.y/sph->world_size.y;
	sim_ratio.z=real_world_side.z/sph->world_size.z;
}/*
void drawSphere(double r, int lats, int longs) {
     int i, j;
	 float M_PI=3.14;
      for(i = 0; i <= lats; i++) {
           double lat0 = 3.14 * (-0.5 + (double) (i - 1) / lats);
           double z0  = sin(lat0);
          double zr0 =  cos(lat0);
    
         double lat1 = M_PI * (-0.5 + (double) i / lats);
           double z1 = sin(lat1);
           double zr1 = cos(lat1);
    
         //  glBegin(GL_QUAD_STRIP);
           for(j = 0; j <= longs; j++) {
               double lng = 2 * M_PI * (double) (j - 1) / longs;
               double x = cos(lng);
               double y = sin(lng);
    
               glNormal3f(x * zr0, y * zr0, z0);
               glVertex3f(x * zr0, y * zr0, z0);
               glNormal3f(x * zr1, y * zr1, z1);
               glVertex3f(x * zr1, y * zr1, z1);
			   glScalef ( 0.5, 0.5, 0.5 );
           }
          // glEnd();
       }
   }*/
void drawSphere(GLfloat xx, GLfloat yy, GLfloat zz, GLfloat radius, GLfloat M, GLfloat N)  
{  

 float step_z = PI/M;  
 float step_xy = 2*PI/N;  
 float x[4],y[4],z[4];  
  
 float angle_z = 0.0;  
 float angle_xy = 0.0;  
 int i=0, j=0;  
 
 glBegin(GL_QUADS);
  for(i=0; i<M; i++)  
  {  
   angle_z = i * step_z;  
     
   for(j=0; j<N; j++)  
   {  
    angle_xy = j * step_xy;  
  
    x[0] = radius * sin(angle_z) * cos(angle_xy);  
    y[0] = radius * sin(angle_z) * sin(angle_xy);  
    z[0] = radius * cos(angle_z);  
  
    x[1] = radius * sin(angle_z + step_z) * cos(angle_xy);  
    y[1] = radius * sin(angle_z + step_z) * sin(angle_xy);  
    z[1] = radius * cos(angle_z + step_z);  
  
    x[2] = radius*sin(angle_z + step_z)*cos(angle_xy + step_xy);  
    y[2] = radius*sin(angle_z + step_z)*sin(angle_xy + step_xy);  
    z[2] = radius*cos(angle_z + step_z);  
  
    x[3] = radius * sin(angle_z) * cos(angle_xy + step_xy);  
    y[3] = radius * sin(angle_z) * sin(angle_xy + step_xy);  
    z[3] = radius * cos(angle_z);  
  
    for(int k=0; k<4; k++)  
    {  
     glVertex3f(xx+x[k], yy+y[k],zz+z[k]);  
    }  
   }  
  }  

 glEnd();  
}  
void render_particles()
{
	glColor3f(0.0f, 0.0f, 1.0f);
	for(uint i=0; i<sph->num_particle; i++)
	{
		sph->mem[i].CalcParticleColor();
		    glColor3f(sph->mem[i].particle_color.x, sph->mem[i].particle_color.y, sph->mem[i].particle_color.z);
			drawSphere(sph->mem[i].pos.x*sim_ratio.x+real_world_origin.x, 
				sph->mem[i].pos.y*sim_ratio.y+real_world_origin.y,
				sph->mem[i].pos.z*sim_ratio.z+real_world_origin.z, 
				0.3,12, 12)  ;	
	}
}

void display_func()
{
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	/////////////////////////////back_ground_color/////////////////////
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	/////////////////////////////////
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glShadeModel(GL_SMOOTH);
    glPushMatrix();

	if(buttonState == 1)
	{
		xRot+=(xRotLength-xRot)*0.1f;
		yRot+=(yRotLength-yRot)*0.1f;
	}

	glTranslatef(xTrans, yTrans, zTrans);
    glRotatef(xRot, 1.0f, 0.0f, 0.0f);
    glRotatef(yRot, 0.0f, 1.0f, 0.0f);

	sph->animation();

	glUseProgram(p);
	render_particles();

	glUseProgram(0);
	draw_box(real_world_origin.x, real_world_origin.y, real_world_origin.z, real_world_side.x, real_world_side.y, real_world_side.z);

	glPopMatrix();

    glutSwapBuffers();
	
	sph_timer->update();
	memset(window_title, 0, 50);
	sprintf(window_title, "SPH System 3D. FPS: %f", sph_timer->get_fps());
	glutSetWindowTitle(window_title);
}

void idle_func()
{
	glutPostRedisplay();
}

void reshape_func(GLint width, GLint height)
{
	window_width=width;
	window_height=height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0, (float)width/height, 0.001, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -3.0f);
}

void keyboard_func(unsigned char key, int x, int y)
{
	if(key == ' ')
	{
		sph->sys_running=1-sph->sys_running;
	}

	if(key == 'w')
	{
		zTrans += 0.3f;
	}

	if(key == 's')
	{
		zTrans -= 0.3f;
	}

	if(key == 'a')
	{
		xTrans -= 0.3f;
	}

	if(key == 'd')
	{
		xTrans += 0.3f;
	}

	if(key == 'q')
	{
		yTrans -= 0.3f;
	}

	if(key == 'e')
	{
		yTrans += 0.3f;
	}

	glutPostRedisplay();
}

void mouse_func(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
        buttonState = 1;
	}
    else if (state == GLUT_UP)
	{
        buttonState = 0;
	}

    ox = x; oy = y;

    glutPostRedisplay();
}

void motion_func(int x, int y)
{
    float dx, dy;
    dx = (float)(x - ox);
    dy = (float)(y - oy);

	if (buttonState == 1) 
	{
		xRotLength += dy / 5.0f;
		yRotLength += dx / 5.0f;
	}

	ox = x; oy = y;

	glutPostRedisplay();
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE| GLUT_DEPTH);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("SPH Fluid 3D");

	init_sph_system();
	init();
	init_ratio();
	//set_shaders();
	
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
	glEnable(GL_POINT_SPRITE_ARB);
	glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(display_func);
	glutReshapeFunc(reshape_func);
	glutKeyboardFunc(keyboard_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutIdleFunc(idle_func);

    glutMainLoop();

    return 0;
}
