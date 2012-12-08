/*
CSCI 480
Assignment 3 Raytracer

Name: Shengzhe Chen
*/

#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>

// define the boundaries of certain variables
#define MAX_TRIANGLES 10
#define MAX_SPHERES 10
#define MAX_LIGHTS 10
#define RANDOM_LIGHTS 20
#define LIGHT_VOLUME 0.05f
#define SAMPLING_TIMES 9


//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];
char *filename=0;

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
  double panel[4];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void calculatePanelEQ(double p1[3], double p2[3], double p3[3], double EQ[4]);
void rayTrace(double from[3], double to[3], int color[3], int depth);
double degreeToRadian(double degree);
double normalize(double point[3]);
double dotMultiple(double vect1[3], double vect2[3]);
void calculatePanelEQ(double p1[3], double p2[3], double p3[3], double EQ[4]);
double distance(double p1[3], double p2[3]);
double triangleArea(double v1[3], double v2[3], double v3[3]);
bool isPointInsideTriangle(double point[3], int triIndex);
bool isIntersectedWithTriangle(double from[3], double to[3], double target[3], int index, double *thres);
bool isIntersectedWithSphere(double from[3], double to[3], double target[3], int index, double *thres);
bool isBlockedByObjects(double from[3], double to[3]);
void ReflectVector(double L[3], double N[3], double R[3]);
void computeTriNormal(double point[3], int triIndex, double N[3]);
void computeTriDif(double point[3], int triIndex, double dif[3]);
void computeTriSpec(double point[3], int triIndex, double spec[3]);
void computeTriShin(double point[3], int triIndex, double *shinn);
void phongShading(double target[3], double V[3], double L[3], double N[3], double R[3], bool isSphere, int lIndex, int index, int color[3], int depth, double lightFactor);
bool isIntersected(double from[3], double to[3], double target[3], bool *isSphere, int *index);
bool pointsInLine(double v1[3], double v2[3], double v3[3]);
double randomLight(int lIndex, double from[3]);
void rayTrace(double from[3], double to[3], int color[3], int depth);
void beginTrace();

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
	glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
	glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  char *file = "raytracing.jpg";

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", file);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(file, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(stricmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

//Convert degree of angle to radian value
double degreeToRadian(double degree)
{
	return degree/180.0*3.1415926;
}
       
int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(stricmp(type,"triangle")==0)
	  {
		printf("found triangle\n");
		int j;
		for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  calculatePanelEQ(t.v[0].position, t.v[1].position, t.v[2].position, t.panel);
	  triangles[num_triangles++] = t;
	}
    else if(stricmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	  {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	  }
	  spheres[num_spheres++] = s;
	}
    else if(stricmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
    else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
  }
  return 0;
}

// normalize a vector
double normalize(double point[3])
{
	double length = sqrt((pow(point[0], 2)+pow(point[1], 2)+pow(point[2], 2)));
	for (int i=0;i<3;i++)
	{
		point[i]/=length;
	}
	return length;
}

// dot product of two vectors
double dotMultiple(double vect1[3], double vect2[3])
{
	double ret = 0;
	for (int i=0;i<3;i++)
	{
		ret += vect1[i]*vect2[i];
	}
	return ret;
}

// calculate the panel equation based on the points of triangle
void calculatePanelEQ(double p1[3], double p2[3], double p3[3], double EQ[4])
{
	// calculate the normal of the panel and then compute the last cof of the panel equation

	double vector1[3], vector2[3], normal[3];
	for (int i=0;i<3;i++)
	{
		vector1[i] = p2[i]-p1[i];
		vector2[i] = p2[i]-p3[i];
	}
	normal[0] = vector1[1]*vector2[2]-vector1[2]*vector2[1];
	normal[1] = vector1[2]*vector2[0]-vector1[0]*vector2[2];
	normal[2] = vector1[0]*vector2[1]-vector1[1]*vector2[0];
	normalize(normal);
	
	for (int i=0;i<3;i++)
	{
		EQ[i] = normal[i];
	}

	EQ[3] = -1*(normal[0]*p1[0]+normal[1]*p1[1]+normal[2]*p1[2]);
}

// calculate the distance between two points
double distance(double p1[3], double p2[3])
{
	double ret = 0.0;
	for (int i=0;i<3;i++)
	{
		ret += pow((p1[i]-p2[i]), 2);
	}

	return sqrt(ret);
}

// calculate the area of a triangle
double triangleArea(double v1[3], double v2[3], double v3[3])
{
	// compute the triangle area using S= 1/4*sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a))
	double area = 0.0;

	double a = distance(v1, v2);
	double b = distance(v1, v3);
	double c = distance(v2, v3);
	double delta = (a+b+c)*(a+b-c)*(a+c-b)*(b+c-a);
	if (delta>0)
	{
		area = 0.25*sqrt((a+b+c)*(a+b-c)*(a+c-b)*(b+c-a));
	}
	
	return area;
}

// check whether a point is inside a triangle 
bool isPointInsideTriangle(double point[3], int triIndex)
{
	// if point is inside the triangle, the area of the triangles that has the 4th point should be equal to the original triangle

	Vertex v1, v2, v3;
	v1 = triangles[triIndex].v[0];
	v2 = triangles[triIndex].v[1];
	v3 = triangles[triIndex].v[2];

	double v12[3], v13[3], triNormal[3];
	for (int i=0;i<3;i++)
	{
		v12[i] = v2.position[i]-v1.position[i];
		v13[i] = v3.position[i]-v1.position[i];
	}
	triNormal[0] = v12[1]*v13[2]-v12[2]*v13[1];
	triNormal[1] = v12[2]*v13[0]-v12[0]*v13[2];
	triNormal[2] = v12[0]*v13[1]-v12[1]*v13[0];
	normalize(triNormal);
	if (triNormal[2] != 0)
	{
		v1.position[2] = v2.position[2] = v3.position[2] = point[2] = 0;
	}
	else 
	{
		v1.position[1] = v2.position[1] = v3.position[1] = point[1] = 0;
	}
	double triArea = triangleArea(v1.position, v2.position, v3.position);
	double triArea1 = triangleArea(point, v1.position, v2.position);
	double triArea2 = triangleArea(point, v1.position, v3.position);
	double triArea3 = triangleArea(point, v2.position, v3.position);
	double area = triArea1+triArea2+triArea3;

	if (area <= triArea*1.01 && area >= triArea*0.95)
	{
		return true;
	}
	return false;
}

// check whether a ray hits triangle on the space
bool isIntersectedWithTriangle(double from[3], double to[3], double target[3], int index, double *thres)
{
	double d[3];
	for (int i=0;i<3;i++)
	{
		d[i] = to[i]-from[i];
	}
	normalize(d);
	double nd = dotMultiple(triangles[index].panel, d);
	if (nd == 0)
		return false;

	double t = -1*(dotMultiple(triangles[index].panel, from)+triangles[index].panel[3])/nd;
	if (t<=0.01)
		return false;
	
	if (*thres > t || *thres == 0.0)
	{
		double tempTar[3];
		for (int i=0;i<3;i++)
		{
			tempTar[i] = from[i]+t*d[i];
		}


		if (isPointInsideTriangle(tempTar, index))
		{
			for (int i=0; i<3; i++)
			{
				target[i] = from[i]+t*d[i];
			}
			*thres = t;
			return true;
		}
	}
	return false;
}

// check whether a ray hits sphere on the space
bool isIntersectedWithSphere(double from[3], double to[3], double target[3], int index, double *thres)
{
	double d[3];
	for (int i=0;i<3;i++)
	{
		d[i] = to[i]-from[i];
	}
	normalize(d);
	double b = 2*(d[0]*(from[0]-spheres[index].position[0])+d[1]*(from[1]-spheres[index].position[1])+d[2]*(from[2]-spheres[index].position[2]));
	double c = pow((from[0]-spheres[index].position[0]), 2)+pow((from[1]-spheres[index].position[1]), 2)+pow((from[2]-spheres[index].position[2]), 2)-pow(spheres[index].radius, 2);
	double b24c = pow(b, 2)-4*c;
	if (b24c<0.2)
		return false;
		
	double t0 = (-b+sqrt(b24c))/2;
	double t1 = (-b-sqrt(b24c))/2;
	double t = 0.0;

	if (t0<=0 && t1<=0)
		return false;

	if (t1<=0 || t0<=0)
		t = max(t0, t1);
	else 
		t = min(t0, t1);
	if (*thres<t && *thres != 0.0)
		return false;
	
	if (t<0.01)
		return false;
	
	for (int i=0; i<3; i++)
		target[i] = from[i] + t*d[i];
	*thres = t;
	return true;
}

// check whether a ray is blocked by obstacles from light
bool isBlockedByObjects(double from[3], double to[3])
{
	// checking block issue starts from checking spheres and then triangles

	double target[3] = {0, 0, 0};
	double thres = 0.0;
	for (int i=0; i<num_spheres; i++)
	{
		if(isIntersectedWithSphere(from, to, target, i, &thres))
		{
			if (distance(from, target) <=distance(from, to))
			{
				return true;
			}
		}
	}

	for (int i=0; i<num_triangles; i++)
	{
		if (isIntersectedWithTriangle(from, to, target, i, &thres))
		{
			if (distance(from, target) <= distance(from, to))
			{
				return true;
			}
		}
	}

	return false;
}

// calculate the reflected vector based on the view, and normal vector
void ReflectVector(double L[3], double N[3], double R[3])
{
	double ln = dotMultiple(L, N);
	if (ln<0)
		ln = 0;
	 
	for (int i=0; i<3; i++)
	{
		R[i] = 2*ln*N[i] - L[i];
	}
	normalize(R);
}

// interpolate the normal of a specific point inside a triangle
void computeTriNormal(double point[3], int triIndex, double N[3])
{
	double area = triangleArea(triangles[triIndex].v[0].position, triangles[triIndex].v[1].position, triangles[triIndex].v[2].position);
	double area01 = triangleArea(point, triangles[triIndex].v[0].position, triangles[triIndex].v[1].position);
	double area02 = triangleArea(point, triangles[triIndex].v[0].position, triangles[triIndex].v[2].position);
	double area12 = triangleArea(point, triangles[triIndex].v[1].position, triangles[triIndex].v[2].position);
	double f01 = area01/area;
	double f02 = area02/area;
	double f12 = area12/area;

	for (int i=0; i<3; i++)
	{
		N[i] = f01*triangles[triIndex].v[2].normal[i]+f02*triangles[triIndex].v[1].normal[i]+f12*triangles[triIndex].v[0].normal[i];
	}

	normalize(N);
}

// interpolate the diffuse cof of a specific point inside a triangle
void computeTriDif(double point[3], int triIndex, double dif[3])
{
	double area = triangleArea(triangles[triIndex].v[0].position, triangles[triIndex].v[1].position, triangles[triIndex].v[2].position);
	double area01 = triangleArea(point, triangles[triIndex].v[0].position, triangles[triIndex].v[1].position);
	double area02 = triangleArea(point, triangles[triIndex].v[0].position, triangles[triIndex].v[2].position);
	double area12 = triangleArea(point, triangles[triIndex].v[1].position, triangles[triIndex].v[2].position);
	double f01 = area01/area;
	double f02 = area02/area;
	double f12 = area12/area;

	for (int i=0; i<3; i++)
	{
		dif[i] = f01*triangles[triIndex].v[2].color_diffuse[i]+f02*triangles[triIndex].v[1].color_diffuse[i]+f12*triangles[triIndex].v[0].color_diffuse[i];
		if (dif[i]>1)
			dif[i] = 1.0;
		else if (dif[i]<0)
			dif[i] = 0;
	}
}

// interpolate the specular cof of a specific point inside a triangle
void computeTriSpec(double point[3], int triIndex, double spec[3])
{
	double area = triangleArea(triangles[triIndex].v[0].position, triangles[triIndex].v[1].position, triangles[triIndex].v[2].position);
	double area01 = triangleArea(point, triangles[triIndex].v[0].position, triangles[triIndex].v[1].position);
	double area02 = triangleArea(point, triangles[triIndex].v[0].position, triangles[triIndex].v[2].position);
	double area12 = triangleArea(point, triangles[triIndex].v[1].position, triangles[triIndex].v[2].position);
	double f01 = area01/area;
	double f02 = area02/area;
	double f12 = area12/area;

	for (int i=0; i<3; i++)
	{
		spec[i] = f01*triangles[triIndex].v[2].color_specular[i]+f02*triangles[triIndex].v[1].color_specular[i]+f12*triangles[triIndex].v[0].color_specular[i];
		if (spec[i]>1)
			spec[i] = 1.0;
		else if (spec[i]<0)
			spec[i] = 0;
	}
}

// interpolate the shinn of a specific point inside a triangle
void computeTriShin(double point[3], int triIndex, double *shinn)
{
	double area = triangleArea(triangles[triIndex].v[0].position, triangles[triIndex].v[1].position, triangles[triIndex].v[2].position);
	double area01 = triangleArea(point, triangles[triIndex].v[0].position, triangles[triIndex].v[1].position);
	double area02 = triangleArea(point, triangles[triIndex].v[0].position, triangles[triIndex].v[2].position);
	double area12 = triangleArea(point, triangles[triIndex].v[1].position, triangles[triIndex].v[2].position);
	double f01 = area01/area;
	double f02 = area02/area;
	double f12 = area12/area;
	*shinn = f01*triangles[triIndex].v[2].shininess+f02*triangles[triIndex].v[1].shininess+f12*triangles[triIndex].v[0].shininess;

	if (*shinn<0)
		*shinn = 0;
}

// local phong shading
void phongShading(double target[3], double V[3], double L[3], double N[3], double R[3], bool isSphere, int lIndex, int index, int color[3], int depth, double lightFactor)
{
	// different intersected objects should have different shading mechnism based on their normals and other details

	if (isSphere)
	{
		double ln = dotMultiple(L, N);
		if (ln<0)
			ln = 0;
		double rv = dotMultiple(R, V);
		if (rv<0)
			rv = 0;
		double rva = pow(rv, spheres[index].shininess);
		for (int i=0; i<3; i++)
		{
			color[i] += (int)((ln*spheres[index].color_diffuse[i]+rva*spheres[index].color_specular[i])*lights[lIndex].color[i]*255*lightFactor);
		}
	}
	else 
	{
		double ln = dotMultiple(L, N);
		if (ln<0)
			ln = 0;
		double rv = dotMultiple(R, V);
		if (rv<0)
			rv = 0;
		double spec[3], dif[3], shinn;
		computeTriDif(target, index, dif);
		computeTriSpec(target, index, spec);
		computeTriShin(target, index, &shinn);
		double rva = pow(rv, shinn);
		for (int i=0; i<3; i++)
		{
			color[i] += (int)((ln*dif[i]+rva*spec[i])*lights[lIndex].color[i]*255*lightFactor);
		}
	}
}

// check whether a ray is intersected with any objects on the scene
bool isIntersected(double from[3], double to[3], double target[3], bool *isSphere, int *index)
{
	bool isIntersected = false;
	double thres = 0.0;
	for (int i=0; i<num_spheres; i++)
	{
		if(isIntersectedWithSphere(from, to, target, i, &thres))
		{
			isIntersected = true;
			*index = i;
			*isSphere = true;
		}
	}
	for (int i=0; i<num_triangles; i++)
	{
		if(isIntersectedWithTriangle(from, to, target, i, &thres))
		{
			isIntersected = true;
			*index = i;
			*isSphere = false;
		}
	}

	return isIntersected;
}

// check whether three arbitrary points are on the same line
bool pointsInLine(double v1[3], double v2[3], double v3[3])
{
	double v12[3], v13[3];
	for (int i=0; i<3; i++)
	{
		v12[i] = v2[i]-v1[i];
		v13[i] = v3[i]-v1[i];
	}
	normalize(v12);
	normalize(v13);

	for (int i=0; i<3; i++)
	{
		if (fabs(v12[i]-v13[i]) > 0.00001)
		{
			return false;
		}
	}
	return true;
}

// point light is regarded as a volume light source, this function randomly computes light contribution factor in order to have soft shadow effect
double randomLight(int lIndex, double from[3])
{
	// randomly pick up a point surrounding the light center as the light position to get how much the light affects the visual effect
	double lightPosition[3] = {0, 0, 0};
	int count = 0;
	for (int i=0; i<RANDOM_LIGHTS; i++)
	{
		double degree1 = (rand()%360)*1.0;
		double degree2 = (rand()%360)*1.0;
		lightPosition[0] = LIGHT_VOLUME*cos(degreeToRadian(degree1))*sin(degreeToRadian(degree2)) + lights[lIndex].position[0];
		lightPosition[1] = LIGHT_VOLUME*sin(degreeToRadian(degree1))*sin(degreeToRadian(degree2)) + lights[lIndex].position[1];
		lightPosition[2] = LIGHT_VOLUME*cos(degreeToRadian(degree2)) + lights[lIndex].position[2];

		if (!isBlockedByObjects(from, lightPosition))
			count++;
	}

	return count*1.0/RANDOM_LIGHTS;
}

// follow a ray shotting from the eye
void rayTrace(double from[3], double to[3], int color[3], int depth)
{
	// a ray is limited to be kept tracing for more than 5 times
	if (depth>=5)
		return;

	int index = -1;
	bool isSphere = false;
	double target[3] = {0, 0, 0};

	bool isIntered = isIntersected(from, to, target, &isSphere, &index);
	double V[3], L[3], N[3], R[3];

	for (int i=0;i<num_lights;i++)
	{
		if (pointsInLine(from, to, lights[i].position) && (distance(from, lights[i].position)<distance(target, lights[i].position)))
		{
			for (int j=0; j<3; j++)
			{
				color[j] = lights[i].color[j];
			}
			return;
		}
	}

	if (isIntered)
	{
		for (int i=0; i<num_lights; i++)
		{
			double lightFactor = randomLight(i, target);
			if (lightFactor>0.3)
			{
				if (isSphere)
				{
					for (int j=0; j<3; j++)
					{
						V[j] = from[j] - target[j];
						L[j] = lights[i]. position[j]- target[j];
						N[j] = target[j] - spheres[index].position[j];
					}
					normalize(V);
					normalize(L);
					normalize(N);
					ReflectVector(L, N, R);
					phongShading(target, V, L, N, R, isSphere, i, index, color, depth, lightFactor);
					
					int reflectionColor[3] = {0, 0, 0};
					double tempTo[3] = {0, 0, 0};
					for (int j=0; j<3; j++)
					{
						tempTo[j] = target[j] + R[j];
					}
					rayTrace(target, tempTo, reflectionColor, depth+1);
					for (int j=0; j<3; j++)
					{
						color[j] += reflectionColor[j]*pow(spheres[index].color_specular[j], depth+1);
					}
				}
				else
				{
					for (int j=0; j<3; j++)
					{
						V[j] = from[j] - target[j];
						L[j] = lights[i]. position[j]- target[j];
					}
					normalize(V);
					normalize(L);
					computeTriNormal(target, index, N);
					ReflectVector(L, N, R);
					phongShading(target, V, L, N, R, isSphere, i, index, color, depth, lightFactor);
					
					int reflectionColor[3] = {0, 0, 0};
					double tempTo[3] = {0, 0, 0};
					for (int j=0; j<3; j++)
					{
						tempTo[j] = target[j] + R[j];
					}
					rayTrace(target, tempTo, reflectionColor, depth+1);
					double spec[3] = {0, 0, 0};
					computeTriSpec(target, index, spec);
					for (int j=0; j<3; j++)
					{
						color[j] += reflectionColor[j]*pow(spec[j], depth+1);
					}
				}
			}
			else 
			{
				if (isSphere)
				{
					for (int j=0; j<3; j++)
					{
						V[j] = from[j] - target[j];
						L[j] = lights[i]. position[j]- target[j];
						N[j] = target[j] - spheres[index].position[j];
					}
					normalize(V);
					normalize(L);
					normalize(N);
					ReflectVector(L, N, R);
					int reflectionColor[3] = {0, 0, 0};
					double tempTo[3] = {0, 0, 0};
					for (int j=0; j<3; j++)
					{
						tempTo[j] = target[j] + R[j];
					}
					rayTrace(target, tempTo, reflectionColor, depth+1);
					for (int j=0; j<3; j++)
					{
						color[j] += reflectionColor[j]*pow(spheres[index].color_specular[j], depth+1);
					}
				}
				else
				{
					for (int j=0; j<3; j++)
					{
						V[j] = from[j] - target[j];
						L[j] = lights[i]. position[j]- target[j];
					}
					normalize(V);
					normalize(L);
					computeTriNormal(target, index, N);
					ReflectVector(L, N, R);
					int reflectionColor[3] = {0, 0, 0};
					double tempTo[3] = {0, 0, 0};
					for (int j=0; j<3; j++)
					{
						tempTo[j] = target[j] + R[j];
					}
					rayTrace(target, tempTo, reflectionColor, depth+1);
					double spec[3] = {0, 0, 0};
					computeTriSpec(target, index, spec);
					for (int j=0; j<3; j++)
					{
						color[j] += reflectionColor[j]*pow(spec[j], depth+1);
					}
				}
			}
		}
		for (int i=0; i<3; i++)
		{
			if (depth == 0)
			{
				color[i] += (int)(ambient_light[i]*256);
			}
		}
	}
	else
	{
		color[0] = color[1] = color[2] = 255;
	}

	for (int i=0; i<3;i++)
	{
		if (color[i]>255)
			color[i] = 255;
		else if(color[i]<0)
			color[i] = 0;
	}

	return;
}

// ray tracing main route
void beginTrace()
{
	// in order to have better result, using anti-alias mechnism here, simply choose a given number sampling points around each pixel and then value-average the color to get the final color which should be sent to color buffer
	double aspect = 1.0*WIDTH/HEIGHT;
	double xmin = -aspect*tan(degreeToRadian(fov/2.0));
	double ymin = -tan(degreeToRadian(fov/2.0));
	double xstep = fabs(2.0*xmin)/WIDTH;
	double ystep = fabs(2.0*ymin)/HEIGHT;

	double from[3]={0.0, 0.0, 0.0};
	double to[3];
	to[2]=-1;
	int i, j;
	for (double y=ymin, j=0;y<=fabs(ymin);y+=ystep,j++)
	{
		glBegin(GL_POINTS);
		for (double x=xmin, i=0;x<=fabs(xmin);x+=xstep, i++)
		{
			to[0]=x;
			to[1]=y;
			int finalColor[3] = {0, 0, 0};
			for (int times = 0; times<SAMPLING_TIMES; times++)
			{
				to[0] = x-xstep*0.5f+(rand()%11)*0.1*xstep;
				to[1] = y-ystep*0.5+(rand()%11)*0.1*ystep;
				int color[3] = {0, 0, 0};
				rayTrace(from, to, color, 0);
				for (int fi=0;fi<3;fi++)
				{
					finalColor[fi]+=color[fi];
				}
			}
			for (int ti=0;ti<3;ti++)
			{
				finalColor[ti] = (int)(1.0*finalColor[ti]/SAMPLING_TIMES);
			}
			plot_pixel(i,j,finalColor[0], finalColor[1], finalColor[2]);
		}	
		glEnd();
		glFlush();

	}
}

void init()
{
	glMatrixMode(GL_PROJECTION);
	glOrtho(0,WIDTH,0,HEIGHT,1,-1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{  
	static int once=0;
	if(!once)
	{
		beginTrace();
		if(mode == MODE_JPEG)
			save_jpg();
	}
	once=1;
}

void display()
{
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
   {
	   mode = MODE_JPEG;
	   filename = argv[2];
	}
  else if(argc == 2)
    mode = MODE_DISPLAY;
	
  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
