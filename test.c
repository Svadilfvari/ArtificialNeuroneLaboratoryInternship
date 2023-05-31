#include "myBmpGris.h"

#define IMAGE_NAME "Test"


/*int** greyMatrice(BmpImg img)       //peut-être pas nécessaire
{
    int** matrice[img->dimX][img->dimY];        //allocation dynamique?
    for(int i=0;i<img->dimX;i++)
    {
        for(int j=0;j<img->dimX;j++)
        {
            matrice[i][j]=getPixel(img,i,j);
        }
    }
    return matrice;
}*/

int** diffMatrice(BmpImg img1, BmpImg img2)
{
    if (img1->dimX!=img2->dimX || img1->dimY!=img2->dimY) {return NULL;}
    int** matrice = calloc(img1->dimX,sizeof(int*));
    for(int i=0;i<img1->dimX;i++)
    {
        matrice[i] = calloc(img1->dimY,sizeof(int));
        for(int j=0;j<img1->dimY;j++)
        {
            matrice[i][j]=mat2[i][j]-mat1[i][j];
        }
    }
    return matrice;
}


int main()
{
    BmpImg img1 = readBmpImage(IMAGE_NAME+"1.bmp");
    BmpImg img2 = readBmpImage(IMAGE_NAME+"2.bmp");
    
    int** diffMatrice = diffMatrice(img1,img2);


    if (diffMatrice!=NULL)
    {
        for(int i=0;i<img1->dimX;i++) {free(diffMatrice[i]);}
        free(diffMatrice);
    } diffMatrice=NULL;
    freeBmpImg(img1); img1=NULL;
    freeBmpImg(img2); img2=NULL;

    return 0;
}