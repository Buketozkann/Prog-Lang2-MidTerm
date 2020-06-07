#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct{
	double x;
	double y;
	double z;
} Vector;

void print_vector(const Vector v){
	printf("(%.2lf, %.2lf, %.2lf)",v.x,v.y,v.z);
};
Vector sum(const Vector v1, const Vector v2) {	
	Vector v;
	
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	v.z = v1.z + v2.z;
	
	return v;
};
Vector diff(const Vector v1, const Vector v2){
	Vector v;
	
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	v.z = v1.z - v2.z;
	
	return v;
};

Vector multiplyby_scalar(const Vector v1, const double k) {
	Vector v;
	
	v.x = v1.x * k;
	v.y = v1.y * k;
	v.z = v1.z * k;
	
	return v;
};
double dot_product(const Vector v1, const Vector v2) {
	double product = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	return product;

};
Vector cross_product(const Vector v1, const Vector v2) {
	Vector v;
	
	v.x = v1.y*v2.z - v1.z*v2.y;
	v.y = v1.z*v2.x - v1.x*v2.z;
	v.z = v1.x*v2.y - v1.y*v2.x;
	
	return v;
};
double norm(const Vector v) {
	double a = sqrt(pow(v.x,v.x) + pow(v.y,v.y) + pow(v.z,v.z));
	
	return a;
};
int is_unitvector(const Vector v) {
	double b = sqrt(pow(v.x,v.x) + pow(v.y,v.y) + pow(v.z,v.z));
	
	if (b == 1) {
		return 1;
	}
	else {
		return 0;
	}
};
Vector unit(const Vector v) {
	Vector s; 
	int	r = sqrt(pow(v.x,v.x) + pow(v.y,v.y) + pow(v.z,v.z));
	s.x = v.x/r;
	s.y = v.y/r;
	s.z = v.z/r;
	return s;
};

double angle(const Vector v1, const Vector v2) {
	
  double normone=sqrt(pow(v1.x,2)+pow(v1.y,2)+pow(v1.z,2));
  double normtwo=sqrt(pow(v2.x,2)+pow(v2.y,2)+pow(v2.z,2));
  double dotproduct=v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  double angle = acos(dotproduct/(normone*normtwo));
	return angle;
};

double distance(const Vector v1, const Vector v2) {
	Vector v;
	
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	v.z = v1.z - v2.z;
	
	double dist = sqrt(pow(v.x,2)+pow(v.y,2)+pow(v.z,2));
	
	return dist;
};
int are_linearly_independent(const Vector v1, const Vector v2, const Vector v3) {
	
   int determinant = v1.x * ((v2.y *v3.z) - (v3.y*v1.z -v1.y)) * (v2.x * v3.z - v3.x * v2.z) + v1.z * (v2.x * v3.y - v3.x * v2.y);
	if(determinant == 0)
		return 0;
	else
		return 1;
};
int are_orthogonal(const Vector v1, const Vector v2, const Vector v3) {
	int v1_2 = v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
	int v1_3 = v1.x*v3.x + v1.y*v3.y + v1.z*v3.z;
	int v2_3 = v2.x*v3.x + v2.y*v3.y + v2.z*v3.z;
	
	if ( v1_2 == 0 && v1_3 == 0 && v2_3 == 0)
		return 1;
	else 
		return 0;	
};

int are_orthonormal(const Vector v1,const Vector v2,const Vector v3)
{
  int a = sqrt(pow(v1.x,2)+pow(v1.y,2)+pow(v1.z,2));   
  int b = sqrt(pow(v2.x,2)+pow(v2.y,2)+pow(v2.z,2)); 
  int c = sqrt(pow(v3.x,2)+pow(v3.y,2)+pow(v3.z,2)); 
  int dotproduct_one=v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
  int dotproduct_two=v2.x*v3.x+v2.y*v3.y+ v2.z*v3.z;
  int dotproduct_three=v1.x*v3.x+v1.y*v3.y+v1.z*v3.z;
   if (dotproduct_one==0 && dotproduct_two==0 && dotproduct_three==0 && b==1 && c==1  && a==1 ){  
    return 1;
  }
   else {
    return 0;
   }
}

Vector projection(const Vector v1,const Vector v2)
{
       Vector v;
    
 double product = (v1.x * v2.x )+ (v1.y * v2.y )+ (v1.z * v2.z);
 double n = pow(v2.x,2) + pow(v2.y,2) + pow(v2.z,2);
 	v.x = v2.x * (product/n);
	v.y = v2.y * (product/n);
	v.z = v2.z * (product/n);
	
  return v;
 
};

Vector orthogonal_projection(const Vector v1,const Vector v2)
{
       Vector v;
    
 double product = (v1.x * v2.x )+ (v1.y * v2.y )+ (v1.z * v2.z);
 double n = pow(v2.x,2) + pow(v2.y,2) + pow(v2.z,2);
 	v.x = v1.x - (v2.x * (product/n));
	v.y = v1.y - (v2.y * (product/n));
	v.z = v1.z - (v2.z * (product/n));
	
  return v;
 
};
int convert_2_orthogonal_basis(Vector *v1, Vector *v2, Vector *v3)
{
       Vector u1;
       Vector u2;
       Vector u3;
    u1=*v1;
    u2=diff(*v2,projection(*v2,u1));
    *v2=u2;
    u3=diff(orthogonal_projection(*v3,u1),projection(*v3,u2));
    *v3=u3;
    return 1;
};

char * vector2str(const Vector v)
{
    char *str = (char*)malloc(500*sizeof(char));
    sprintf(str, "%2.lf, %2.lf, %2.lf", v.x, v.y,v.z);
    return str;
    
}

int main () 
{
	Vector v1 = {1, 2, 2}, v2 = {-1, 0, 2}, v3 = {0, 0, 1};
	double k = 2;

    printf("v1 = ");
    print_vector(v1);
    printf("\nv2 = ");
    print_vector(v2);
    printf("\nv3 = ");
    print_vector(v3);
	 
    printf("\nv1 + v2 = ");
	print_vector(sum(v1, v2));
	
	printf("\nV1 - V2 = ");
	print_vector(diff(v1, v2));
	
	printf("\nk * V1 = ");
	print_vector(multiplyby_scalar(v1, k));

	printf("\nV1 . V2 = %.2lf", dot_product(v1, v2));
	
	printf("\nV1 x V2 = ");
	print_vector(cross_product(v1, v2));
    printf("\n| V1 | = %.2lf", norm(v1));
    
    if(is_unitvector(v1))
		printf("\nV1 is a unit vector.");
	else
		printf("\nV1 is not unit vector.");
		
	printf("\nUnit( V1 ) = ");
	print_vector(unit(v1));	
	
	printf("\nangle(V1, V2) = %.2lf", angle(v1, v2));
	
	printf("\ndistance(V1, V2) = %.2lf", distance(v1, v2));
	
	if(are_linearly_independent(v1, v2, v3))
		printf("\nVectors are linearly independent.");
	else
		printf("\nVectors are not linearly independent.");
	if(are_orthogonal(v1, v2, v3))
		printf("\nVectors are orthogonal.");
	else
		printf("\nVectors are not orthogonal.");
	if(are_orthonormal(v1, v2, v3))
		printf("\nVectors are orthonormal.");
	else 
		printf("\nVectors are not orthonormal.");
	
	printf("\nProjection of v1 onto v2 is = ");
    print_vector(projection(v1, v2));
    
    printf("\nOrthogonal Projection of V1 onto V2 is = ");
    print_vector(orthogonal_projection(v1, v2));
   
    if(convert_2_orthogonal_basis(&v1, &v2, &v3)){
		printf("\nOrthogonalization of vectors:");
		printf("\nv1 = ");
	    print_vector(v1);
	    printf("\nv2 = ");
	    print_vector(v2);
	    printf("\nv3 = ");
	    print_vector(v3);
	}	
    printf("\n");
    puts(vector2str(v1)); 
    return 0;  
}
