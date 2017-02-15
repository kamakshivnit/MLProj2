%load proj2matrices.mat
fileID = fopen('Querylevelnorm.txt');
C = textscan(fileID,'%s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s %s  %s %s %s %s %s %s %s %s %s',69623);
%celldisp(C)
fclose(fileID);
%get column vector 3 to 48 store in matrix row
feature = zeros(69623,46); %69623
target = zeros(69623,1);

for i = 1:69623
    target(i,1)= str2double(C{1}{i});
    for j=3:48
      value = strsplit(C{j}{i},':');
      feature(i,j-2)= str2double(value(2));
    end
    disp(i)
end
%disp(feature)
%disp(target)
%save('proj2.mat', 'feature', 'target')
n=69223
random = randperm(n);

sizeT=round(0.8*n);

sizeV=round(0.9*n);

trainSet= feature(random(1:sizeT),:);
validSet= feature(random(sizeT+1:sizeV),:);
testSet= feature(random(sizeV+1:n),:);

targetTrain= target(random(1:sizeT),:);
targetValid= target(random(sizeT+1:sizeV),:);
targetTest= target(random(sizeV+1:n),:);


trainInd1=random(1:sizeT);
trainInd1 = trainInd1(:);
validInd1=random(sizeT+1:sizeV);
validInd1 = validInd1(:);
 % for close form of real data
 
 %commented this hardcoded part for grid search
%M1=3; %hardcoded

%lambda1=0.1
D=46;
validPer1Arr=zeros(10,9);
trainPer1Arr=zeros(10,9);
vaildPer1Min=Inf;
for M1=3:10
    sigmar=0.1*var(trainSet);

    sigma=zeros(D,D);
    sigma=diag(sigmar);
    sigma=0.001*eye(D,D)+sigma;
    %Sigma= diag(sigma)
    Sigma1=zeros(D,D,M1);
    for i=1:M1
     Sigma1(:,:,i)=sigma;
    end

    mu1=zeros(M1,D);
    
    %mu1(:,:)=trainSet(501:503,:); %hardcoded
    mu1Ind=randperm(sizeT,M1);
    mu1(:,:)=trainSet(mu1Ind(1:M1),:);
    mu1=mu1.';


    DesgMatTr=zeros(sizeT,M1);
    for i=1:sizeT
       DesgMatTr(i,1)=1;
    end
    for i=1:sizeT
     for j=2:M1
        DesgMatTr(i,j)=exp(-0.5*(trainSet(i,:).'-mu1(:,j)).'*inv(Sigma1(:,:,j))*(trainSet(i,:).'-mu1(:,j)));
     end
    end
    
    DesgMatVal=zeros(sizeV-sizeT,M1);
    for i=1:sizeV-sizeT
       DesgMatVal(i,1)=1;
    end
    for i=1:sizeV-sizeT
      for j=2:M1
        DesgMatVal(i,j)=exp(-0.5*(validSet(i,:).'-mu1(:,j)).'*inv(Sigma1(:,:,j))*(validSet(i,:).'-mu1(:,j)));
      end
    end

    count=1;
  for lambda1=0.1:0.1:0.9

      w1=inv(lambda1*eye(M1,M1)+(DesgMatTr.'*DesgMatTr))*DesgMatTr.'*targetTrain;

      E=0.5*(targetTrain-DesgMatTr*w1).'*(targetTrain-DesgMatTr*w1);
      trainPer1=sqrt(2*E/sizeT);


      E1=0.5*(targetValid-DesgMatVal*w1).'*(targetValid-DesgMatVal*w1);
      validPer1=sqrt(2*E1/(sizeV-sizeT));
      if (validPer1<vaildPer1Min)
        lambda1min=lambda1;
        M1min=M1;
        trainPer1min=trainPer1;
        vaildPer1Min=validPer1;
        mu1min=mu1;
        w1min=w1;
        sigma1min=Sigma1;
        DesgMatTr_min=DesgMatTr;
      end
      validPer1Arr(M1,count)=validPer1;
      trainPer1Arr(M1,count)=trainPer1;
      count=count+1;
  end
  
end
    lambda1=lambda1min;
    M1=M1min;
    trainPer1=trainPer1min;
    validPer1=vaildPer1Min;
    mu1=mu1min;
    w1=w1min;
    Sigma1=sigma1min;
    DesgMatTr=DesgMatTr_min;
figure
contour3(validPer1Arr)
figure
%close form of synthetic data

load('synthetic.mat')
nsyn=2000
randomsyn = randperm(nsyn);

sizeTsyn=round(0.8*nsyn);

sizeVsyn=round(0.9*nsyn);
x=x.';
x_train= x(randomsyn(1:sizeTsyn),:);
x_valid= x(randomsyn(sizeTsyn+1:sizeVsyn),:);
x_test= x(randomsyn(sizeVsyn+1:nsyn),:);

t_train= t(randomsyn(1:sizeTsyn),:);
t_valid= t(randomsyn(sizeTsyn+1:sizeVsyn),:);
t_test= t(randomsyn(sizeVsyn+1:nsyn),:);


trainInd2=randomsyn(1:sizeTsyn);
trainInd2 = trainInd2(:);
validInd2=randomsyn(sizeTsyn+1:sizeVsyn);
validInd2 = validInd2(:);


%%%%%%%
%M2=3; %hardcoded
D2=10;
validPer2Arr=zeros(10,9);
trainPer2Arr=zeros(10,9);
vaildPer2Min=Inf;
sigmaS=0.1*var(x_train);
for i=1:D2
  if ((sigmaS(i)>-0.00001)&&(sigmaS(i)<0.00001))
      sigmaS(i)=sigmaS(i)+0.1;
  end
end
for M2=3:10
    sigma2=zeros(D2,D2);
    sigma2=diag(sigmaS);
    %Sigma= diag(sigma)
    Sigma2=zeros(D2,D2,M2);
    for i=1:M2
     Sigma2(:,:,i)=sigma2;
    end

    mu2=zeros(M2,D2);

    mu2Ind=randperm(sizeTsyn,M2);
    mu2(:,:)=x_train(mu2Ind(1:M2),:);
    %mu2(:,:)=x_train(101:103,:); %hardcoded
    mu2=mu2.';

    DesgMatTr_S=zeros(sizeTsyn,M2);
    for i=1:sizeTsyn
       DesgMatTr_S(i,1)=1;
    end
    for i=1:sizeTsyn
     for j=2:M2
        DesgMatTr_S(i,j)=exp(-0.5*(x_train(i,:).'-mu2(:,j)).'*inv(Sigma2(:,:,j))*(x_train(i,:).'-mu2(:,j)));
     end
    end

   DesgMatVal_S=zeros(sizeVsyn-sizeTsyn,M2);
    for i=1:sizeVsyn-sizeTsyn
       DesgMatVal_S(i,1)=1;
    end
    for i=1:sizeVsyn-sizeTsyn
     for j=2:M2
        DesgMatVal_S(i,j)=exp(-0.5*(x_valid(i,:).'-mu2(:,j)).'*inv(Sigma2(:,:,j))*(x_valid(i,:).'-mu2(:,j)));
     end
    end

    %lambda2=0.2
    count=1;
    for lambda2=0.1:0.1:0.9
      w2=inv(lambda2*eye(M2,M2)+(DesgMatTr_S.'*DesgMatTr_S))*DesgMatTr_S.'*t_train;

      E2=0.5*(t_train-DesgMatTr_S*w2).'*(t_train-DesgMatTr_S*w2)
      trainPer2=sqrt(2*E2/sizeTsyn);

   
      E3=0.5*(t_valid-DesgMatVal_S*w2).'*(t_valid-DesgMatVal_S*w2);
      validPer2=sqrt(2*E3/(sizeVsyn-sizeTsyn));
      if (validPer2<vaildPer2Min)
        lambda2min=lambda2;
        M1min=M2;
        trainPer2min=trainPer2;
        vaildPer2Min=validPer2;
        mu2min=mu2;
        w2min=w2;
        sigma2min=Sigma2;
        DesgMatTr_Smin=DesgMatTr_S;
      end
      validPer2Arr(M2,count)=validPer2;
      trainPer2Arr(M2,count)=trainPer2;
      count=count+1;
    end

end
  lambda2=lambda2min;
  M2=M1min;
  trainPer2=trainPer2min;
  validPer2=vaildPer2Min;
  mu2=mu2min;
  w2=w2min;
  Sigma2=sigma2min;
  DesgMatTr_S=DesgMatTr_Smin;
figure
contour3(validPer2Arr)
figure

%GRADIENT DESCENT FOR REAL DATA

w01=zeros(M1,1);
w_New=zeros(M1,1);
E1=sizeT;
w01=100*rand(M1,1);
dw1=zeros(M1,E1);
eta1=zeros(1,E1);
eta1(:,:)=1;
w_New=w01;
%for j=1:10
ERMS_OLD=-Inf;

for i=1:E1
 
 dw1(:,i)=eta1(i)*(((targetTrain(i)-w_New.'*DesgMatTr(i,:).')*DesgMatTr(i,:).')-lambda1*w_New);
 w_New=w_New+dw1(:,i);

end

intial = norm((w1- w01),2)
final = norm((w1 - w_New),2)



%GRADIENT DESCENT FOR SYNTHETIC DATA

w02=zeros(M2,1);
w_New2=zeros(M2,1);
J1=sizeTsyn;
w02=100*rand(M2,1);
dw2=zeros(M2,J1);
eta2=zeros(1,J1);
eta2(:,:)=1;
w_New2=w02;
%for j=1:10
ERMS_OLD=-Inf;

for i=1:J1

 dw2(:,i)=eta2(i)*(((t_train(i)-w_New2.'*DesgMatTr_S(i,:).')*DesgMatTr_S(i,:).')-lambda2*w_New2);
 w_New2=w_New2+dw2(:,i);

end

intial = norm((w2- w02),2)
final = norm((w2 - w_New2),2)

%save('proj2.mat M1 mu1 Sigma1 lambda1 trainInd1 validInd1 w1' , '-mat7-binary')
save('proj2.mat', '-mat7-binary');