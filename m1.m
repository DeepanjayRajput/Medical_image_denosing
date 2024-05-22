clear all;
close all;
clc;
%-------------------------------------------------------------------------%
%                         Left Hand side of Algorithm
%-------------------------------------------------------------------------%
[imagename1 imagepath1]=uigetfile('source_images\.jpg;.bmp;.png;.tif;.tiff;.pgm;*.gif','Please choose the first input image');
image_input1=imread(strcat(imagepath1,imagename1));

figure;imshow(image_input1);title('original image');
tic;
%-------------------------------------------------------------------------%
%                                Adding Noise
%-------------------------------------------------------------------------%
X=imnoise(image_input1,'gaussian',0,0.0015);
figure;imshow(X);title('noise image');
%-------------------------------------------------------------------------%
%                                Using Gaussian filter
%-------------------------------------------------------------------------%
[m,n,~]=size(X);
p=m/2;
q=n/2;
d0=0.2;
for i=1:m
  for j=1:n
    distance(i,j)=sqrt(((i-p)^2+(j-q)^2));
  end
end
 distance=mat2gray(distance);
for i=1:m
  for j=1:n
    Hglp(i,j)=(exp(-(distance(i,j))^2/(2*(d0^2))));
  end
end
%-------------------------------------------------------------------------%
%                                Using bilateral filter
%-------------------------------------------------------------------------%
R1=imbilatfilt(X);
figure;imshow(R1);title('after using bilateral');


%-------------------------------------------------------------------------%
%                                Subtracting Noise-Bilateral
%-------------------------------------------------------------------------%
G1 = double(X) - double(R1);
figure;imshow(G1,[]);title('minu image');

G2 = double(X) - double(F);
figure;imshow(G2,[]);title('minu 2 image');
%-------------------------------------------------------------------------%
%                                Using DWT
%-------------------------------------------------------------------------%
XX=size(G1)
[LL,LH,HL,HH]=dwt2(G1,'db1');
figure(5)
subplot(2,2,1);imshow(LL);title('LL band of image');
subplot(2,2,2);imshow(LH);title('LH band of image');
subplot(2,2,3);imshow(HL);title('HL band of image');
subplot(2,2,4);imshow(HH);title('HH band of image');


XY=size(G2)
[A,B,C,D]=dwt2(G2,'db1');
figure(5)
subplot(2,2,1);imshow(A);title('A band of image');
subplot(2,2,2);imshow(B);title('B band of image');
subplot(2,2,3);imshow(C);title('C band of image');
subplot(2,2,4);imshow(D);title('D band of image');

%-------------------------------------------------------------------------%
%                                Doing Thresholding in LH
%-------------------------------------------------------------------------%
figure(6);
tmp = LH;

[m n]= find(LH<26);
for j = 1: length(m)
    tmp(m(j),n(j))=0;
end

[m n]= find(LH>26 & LH <= 230);
for j = 1: length(m)
    tmp(m(j),n(j))=0.8;
end

[m n]= find(LH>230);
for j = 1: length(m)
    tmp(m(j),n(j))=0;
end
subplot(3,3,1);
LHN = im2bw(tmp,0);
imshow(LHN,[]); title('threshoding');


figure(6.0);
tmp0 = B;

[m n]= find(B<26);
for j = 1: length(m)
    tmp0(m(j),n(j))=0;
end

[m n]= find(B>26 & B <= 230);
for j = 1: length(m)
    tmp0(m(j),n(j))=0.8;
end

[m n]= find(B>230);
for j = 1: length(m)
    tmp0(m(j),n(j))=0;
end
subplot(3,3,1);
BN = im2bw(tmp0,0);
imshow(BN,[]); title('2threshoding');

%-------------------------------------------------------------------------%
%                                Doing Thresholding in HL
%-------------------------------------------------------------------------%
tmp1 = HL;

[m n]= find(HL<26);
for j = 1: length(m)
    tmp1(m(j),n(j))=0;
end

[m n]= find(HL>26 & HL <= 230);
for j = 1: length(m)
    tmp1(m(j),n(j))=0.8;
end

[m n]= find(HL>230);
for j = 1: length(m)
    tmp1(m(j),n(j))=0;
end

subplot(3,3,2);
HLN = im2bw(tmp1,0);
imshow(HLN,[]); title('threshoding');


tmp3 = C;

[m n]= find(C<26);
for j = 1: length(m)
    tmp3(m(j),n(j))=0;
end

[m n]= find(C>26 & C <= 230);
for j = 1: length(m)
    tmp3(m(j),n(j))=0.8;
end

[m n]= find(C>230);
for j = 1: length(m)
    tmp3(m(j),n(j))=0;
end

subplot(3,3,2);
CN = im2bw(tmp3,0);
imshow(CN,[]); title('2threshoding');


%-------------------------------------------------------------------------%
%                                Doing Thresholding in HH
%-------------------------------------------------------------------------%
tmp2 =HH ;

[m n]= find(HH<26);
for j = 1: length(m)
    tmp2(m(j),n(j))=0;
end

[m n]= find(HH>26 & HH <= 230);
for j = 1: length(m)
    tmp2(m(j),n(j))=0.8;
end

[m n]= find(HH>230);
for j = 1: length(m)
    tmp2(m(j),n(j))=0;
end
subplot(3,3,3);
HHN = im2bw(tmp2,0);
imshow(HHN,[]); title('threshoding(Between 26-230)');


tmp4 =D ;

[m n]= find(D<26);
for j = 1: length(m)
    tmp4(m(j),n(j))=0;
end

[m n]= find(D>26 & D <= 230);
for j = 1: length(m)
    tmp4(m(j),n(j))=0.8;
end

[m n]= find(D>230);
for j = 1: length(m)
    tmp4(m(j),n(j))=0;
end
subplot(3,3,3);
DN = im2bw(tmp4,0);
imshow(DN,[]); title('2threshoding(Between 26-230)');
%-------------------------------------------------------------------------%
%                                Using IDWT
%-------------------------------------------------------------------------%
figure(7);
R5 = idwt2(LL,LHN,HLN,HHN,'db1',XX);
imshow(R5,[]); title('after idwt2');


figure(7.0);
R7 = idwt2(A,BN,CN,DN,'db1',XY);
imshow(R7,[]); title('2after idwt2');

%-------------------------------------------------------------------------%
%                   Combining IDWT result and Bilateral result
%-------------------------------------------------------------------------%
R6= double(R1) + double(R5);
figure(8);imshow(R6,[]);title('final output image');

R8= double(f) + double(R7);
figure(8.0);imshow(R8,[]);title('2final output image');

toc;