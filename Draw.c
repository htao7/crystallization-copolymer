#include <stdio.h>
#include <stdlib.h>

#define XLength 32
#define YLength 32
#define ZLength 32
#define CHAIN_NUM 1
#define SEGMENT_NUM 32
#define E 3

typedef struct{
	int x, y, z;
}point;

int CoordTrans(int Coord, int EdgeLength);
void InitialArray(FILE *fp);
void POVRayPara(FILE *fp);
void DrawLine(FILE *fp, point p1, point p2);
void DrawCubeFrame(FILE *fp, point p);
void DrawBond(FILE *fp, point p1, point p2);
void DrawArea(FILE *fp);
int jf(int x1, int y1, int z1, int x2, int y2, int z2);

int xx[CHAIN_NUM + 1][SEGMENT_NUM + 2];
int yy[CHAIN_NUM + 1][SEGMENT_NUM + 2];
int zz[CHAIN_NUM + 1][SEGMENT_NUM + 2];
int ii[XLength + 1][YLength + 1][ZLength + 1];
int jj[XLength + 1][YLength + 1][ZLength + 1];

int main (void)
{
  FILE *fp;
  char filename[255];

//  for (float temp = 0.1; temp <= 9.9; temp += 0.1){
	  sprintf(filename, "crystalnun12state0.ord");
	  fp=fopen(filename, "r");
	  if (!fp){
		  printf("Can't open input file!\n");
		  exit(1);
	  }
	  InitialArray(fp);
	  fclose(fp);

	  sprintf(filename, "Grp10Temp%.3fCyc200.pov");
	  fp=fopen(filename, "w");
	  if (!fp) {
		  printf("Can't open output file!\n");
		  exit(1);
	  }
	  POVRayPara(fp);
	  DrawCubeFrame(fp, (point){ XLength*E, YLength*E, ZLength*E });
	  DrawArea(fp);
	  fclose(fp);
 // }
}

int CoordTrans(int Coord, int EdgeLength)
{
	if (Coord < 0)
		return (Coord + (1 - Coord / EdgeLength)*EdgeLength);
	if (Coord % EdgeLength == 0)
		return EdgeLength;
	return Coord % EdgeLength;
}

void InitialArray(FILE *fp)
{
	int temp_x, temp_y, temp_z;
	int i, j, k;
	/* input the initial state*/
	for (i = 1; i <= XLength; i++)
	for (j = 1; j <= YLength; j++)
	for (k = 1; k <= ZLength; k++){
	  ii[i][j][k] = 0;
	  jj[i][j][k] = -1;
	}

	for (i = 1; i <= CHAIN_NUM; i++)
	for (j = 1; j <= SEGMENT_NUM; j++){
	  fscanf(fp, "%i %i %i\n", &temp_x, &temp_y, &temp_z);

	  temp_x = CoordTrans(temp_x, XLength);
	  temp_y = CoordTrans(temp_y, YLength);
	  temp_z = CoordTrans(temp_z, ZLength);

	  xx[i][j] = temp_x;
	  yy[i][j] = temp_y;
	  zz[i][j] = temp_z;

	  ii[temp_x][temp_y][temp_z] = i;
	  jj[temp_x][temp_y][temp_z] = j;
	}
}

void POVRayPara(FILE *fp)
{
	fprintf(fp, "#include \"colors.inc\"\n");
	fprintf(fp, "background {color White}\n");
	fprintf(fp, "global_settings\n{max_trace_level 100}\n");
	fprintf(fp, "camera { location <50,50,-100> angle 30 sky <0, 0, 1> look_at <50, 100, 100> }\n");
	fprintf(fp, "light_source { <-5000,100,100> color White }\n \n");
	fprintf(fp, "light_source { <100,-1000,100> color White }\n \n");
	fprintf(fp, "light_source { <100,100,2000> color White shadowless }\n \n");
	fprintf(fp, "light_source { <100,100,2000> color White shadowless }\n \n");
}

void DrawLine(FILE *fp, point p1, point p2)
{
	fprintf(fp, "cylinder{<%i,%i,%i>,<%i,%i,%i>, 1.0 open \n", p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
	fprintf(fp, "texture {finish {metallic phong 1} pigment {color Black}} rotate <0,0,360*clock>}\n");
}

void DrawCubeFrame(FILE *fp, point p)
{
	DrawLine(fp, (point){  0,   0,   0 }, (point){ p.x,   0,   0 });
	DrawLine(fp, (point){  0, p.y,   0 }, (point){ p.x, p.y,   0 });
	DrawLine(fp, (point){  0,   0, p.z }, (point){ p.x,   0, p.z });
	DrawLine(fp, (point){  0, p.y, p.z }, (point){ p.x, p.y, p.z });
	DrawLine(fp, (point){  0,   0,   0 }, (point){   0, p.y,   0 });
	DrawLine(fp, (point){p.x,   0,   0 }, (point){ p.x, p.y,   0 });
	DrawLine(fp, (point){  0,   0, p.z }, (point){   0, p.y, p.z });
	DrawLine(fp, (point){p.x,   0, p.z }, (point){ p.x, p.y, p.z });
	DrawLine(fp, (point){  0,   0,   0 }, (point){   0,   0, p.z });
	DrawLine(fp, (point){p.x,   0,   0 }, (point){ p.x,   0, p.z });
	DrawLine(fp, (point){  0, p.y,   0 }, (point){   0, p.y, p.z });
	DrawLine(fp, (point){p.x, p.y,   0 }, (point){ p.x, p.y, p.z });
}

void DrawBond(FILE *fp, point p1, point p2)
{
	fprintf(fp, "cylinder{<%i,%i,%i>,<%i,%i,%i>, %f open \n", p1.x*E, p1.y*E, p1.z*E, p2.x*E, p2.y*E, p2.z*E, 0.33*E);
	fprintf(fp, "texture { finish {metallic phong 1} pigment {color rgb <%f,%f,%f>}} rotate <0,0,360*clock>}\n", 30.0 / 100, 30.0 / 100, 30.0 / 100);
	fprintf(fp, "sphere{<%i,%i,%i>, %f\n", p1.x*E, p1.y*E, p1.z*E, 0.43*E);
	fprintf(fp, "texture { finish {metallic phong 1} pigment {color rgb <0.1,0.1,0.9>}} rotate <0,0,360*clock>}\n");
	fprintf(fp, "sphere{<%i,%i,%i>, %f\n", p2.x*E, p2.y*E, p2.z*E, 0.43*E);
	fprintf(fp, "texture { finish {metallic phong 1} pigment {color rgb <0.1,0.1,0.9>}} rotate <0,0,360*clock>}\n");
}

void DrawArea(FILE *fp)
{
	int i, j;
	for (i = 1; i <= CHAIN_NUM; i++){
		for (j = 1; j < SEGMENT_NUM; j++){
			int x1, y1, z1, x2, y2, z2;
			x1 = xx[i][j], x2 = xx[i][j + 1];
			y1 = yy[i][j], y2 = yy[i][j + 1];
			z1 = zz[i][j], z2 = zz[i][j + 1];
			if (abs(x1 - x2)>2 || abs(y1 - y2) > 2 || abs(z1 - z2) > 2)
				continue;
			if (jf(x1, y1, z1, x2, y2, z2) < 0)
				continue;
			DrawBond(fp, (point){ x1, y1, z1 }, (point){ x2, y2, z2 });
			continue;
		}
	}
}

int jf(int x1, int y1, int z1, int x2, int y2, int z2)
{
	int ii1, ii2, jj1, jj2;
	int xs, ys, zs, xe, ye, ze;
	int u, v, w, c = 0;
	x1 = CoordTrans(x1, XLength);
	y1 = CoordTrans(y1, YLength);
	z1 = CoordTrans(z1, ZLength);
	x2 = CoordTrans(x2, XLength);
	y2 = CoordTrans(y2, YLength);
	z2 = CoordTrans(z2, ZLength);
	for (u = -1; u <= 1; u++)
	for (v = -1; v <= 1; v++)
	for (w = -1; w <= 1; w++)
	{ /*三个for循环，将旧键在空间中进行了平移*/
		if (u == 0 && v == 0 && w == 0)
			continue;
		xs = CoordTrans(x1 + u, XLength), ys = CoordTrans(y1 + v, YLength), zs = CoordTrans(z1 + w, ZLength);
		xe = CoordTrans(x2 + u, XLength), ye = CoordTrans(y2 + v, YLength), ze = CoordTrans(z2 + w, ZLength);
		if (xs == x2 && ys == y2 && zs == z2)
			continue;
		if (xe == x1 && ye == y1 && ze == z1)
			continue; /*防止找到的键与原来的键接在一起?*/
		ii1 = ii[xs][ys][zs];
		ii2 = ii[xe][ye][ze];
		jj1 = jj[xs][ys][zs];
		jj2 = jj[xe][ye][ze];
		if (ii1 == ii2 && abs(jj1 - jj2) == 1)
			c++;
	}
	return c;
}
