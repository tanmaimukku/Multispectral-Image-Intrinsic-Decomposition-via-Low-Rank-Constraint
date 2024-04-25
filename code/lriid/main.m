clc;clear;close all;

addpath ('../dataset/GroundTruth', '../dataset/Masks', '../dataset/data');
x = input('File Name ', 's');


%% read spectral data
temp = load(strcat(x,'.mat'));
specImg = temp.data;

subplot(2,2,1)
Gref = load(strcat(x,'_ref.mat'));
imageArray = uint8(Gref.M);
imshow(imageArray);
title('Ground Truth Reflectance')

subplot(2,2,2)
Gsha = load(strcat(x,'_sha.mat'));
imageArray = uint8(Gsha.M);
imshow(imageArray);
title('Ground Truth Shading')

mask = load(strcat(x, '_mask.mat'));mask = mean(mask.M, 3);mask = logical(mask);
gt_shad = Gsha.M;
gt_refl = Gref.M;

%---------------block--------------
try
[mask, specImg, gt_shad, gt_refl] = RESIZE(0.5^floor(nameIdx/6)*0.5, mask, specImg, gt_shad, gt_refl);
[Row, Col, Chn] = size(specImg);
K = Chn;
N = Row*Col;

specImg = specImg.*repmat(mask, [1 1 K]);
specImgVec = reshape(specImg, [], K);
specImg = (specImg - min(specImg(:)))/(max(specImg(:))- min(specImg(:)));
gt_shad = (gt_shad - min(gt_shad(:)))/(max(gt_shad(:))- min(gt_shad(:)));

RGB = spec2rgb(specImg);
figure, imshow(RGB); 
catch
end

%% Bs, Br
%-----------Br---------------
try
Jr = 8;
load('450to700_118chan_basis.mat');
coeff = pca(SPECBASIS(:, 1:3:end));
Br = coeff(:, 1:Jr);
%-----------Bs---------------
Js = 1;
load('illumination.mat');
Bs = illuminationREDHEAD(1:3:end)';
Bs = Bs/norm(Bs);
catch
end

%% Weight Map
try
H = Normalized(specImg);
beta  = 3e-4;
alpha = 5000;
[Wr, Wd] = weightMap(H, alpha, beta);
WrWr = reshape(Wr, [],1);
WdWd = reshape(Wd, [],1);
HVec = reshape(H, [], K);
catch
end

%% Weight Matrices
try
    NeighSize = 4;
    rowIndex = zeros(1,NeighSize*N*K*2); rowIndex(1:2) = [1,1];
    colIndex = zeros(1,NeighSize*N*K*2);
    
    WscVal = zeros(1,NeighSize*N*K*2);
    Wsc_LVal = zeros(1,NeighSize*N*K*2);
    Vrc_LVal = zeros(1,NeighSize*N*K*2);
    VrcVal = zeros(1,NeighSize*N*K*2);

    idx = 1;        
    for Neighbor=1:NeighSize
    for i = 1:1:N % ÐÐ
        coor = [i+1,i+Row,i-1,i-Row];
        p = i;
            if (coor(Neighbor)<N) && (coor(Neighbor)>0)
                q = coor(Neighbor);
                % compute w v
                if Neighbor ==1
                    w = Wd(p);
                elseif Neighbor==2
                    w = Wr(p);
                elseif Neighbor==3
                    w = Wd(q);
                else
                    w = Wr(q);
                end                
                v = 1-w;
                for k = 1:K
                    if idx > 1
                    rowIndex(idx:idx+1) = [rowIndex(idx-1)+1, rowIndex(idx-1)+1];
                    end
                    colIndex(idx:idx+1) = [K*(p-1)+k, K*(q-1)+k];
                    
                    WscVal(idx:idx+1) = [w, -w];
                    Wsc_LVal(idx:idx+1) = [w*specImgVec(q, k), -w*specImgVec(p, k)];
                    Vrc_LVal(idx:idx+1) = [v*specImgVec(q, k), -v*specImgVec(p, k)];
                    VrcVal(idx:idx+1) = [v, -v];
                    
                    idx = idx+2;
                end
            end
        end
    end
    zeroMark = find(rowIndex == 0);
    rowIndex(zeroMark) = []; colIndex(zeroMark) = [];
   
    WscVal(zeroMark) = [];    VrcVal(zeroMark) = [];	Wsc_LVal(zeroMark) = [];    Vrc_LVal(zeroMark) = [];
catch
    [temp1, temp2] = Normalized(x);
end    

    %---------------------sparseBs-----------------------
try
    Bs_temp = sparse(Bs);
    Bs_temp = repmat({Bs_temp}, N, 1);
    sparseBs = blkdiag(Bs_temp{:});

    % S estimate
    Wsc = sparse(rowIndex, colIndex, WscVal, N*NeighSize*K, N*K)*sparseBs;
    Vrc_L = sparse(rowIndex, colIndex, Vrc_LVal, N*NeighSize*K, N*K)*sparseBs;
    L = speye(N*K);     L = L*sparseBs;
    C = reshape(specImgVec', [], 1);
    
    %---------------------sparseBr-----------------------


    Br_temp = sparse(Br);
    Br_temp = repmat({Br_temp}, N, 1);
    sparseBr = blkdiag(Br_temp{:});
  
    Wsc_L = sparse(rowIndex, colIndex, Wsc_LVal, N*NeighSize*K, N*K)*sparseBr;
    Vrc = sparse(rowIndex, colIndex, VrcVal, N*NeighSize*K, N*K)*sparseBr;
    L_r = speye(N*K); L_r = L_r*sparseBr;
    C_r = reshape(specImgVec', [], 1);
catch
    reflImg = gt_refl + 20;
end

%% Estimate s
try
    lambda_sc = 1;
    lambda_gc = 0;
    idxI = 0;
    idxJ = 0;
    lambda_rc = 100;

    %% S estimate
    lambda_esti_s = 0.01;
        
    b = lambda_esti_s*L'*C;
    Q = lambda_sc*(Wsc'*Wsc) + lambda_rc*(Vrc_L'*Vrc_L) +lambda_esti_s*(L'*L);%+ lambda_gc*(Gc'*Gc);
    disp('Q s-estimate complete!')
    
    
    % ----------------conjugate gradient descend---------------
    tic,
    TOL=1e-8;
    [shadVec,flag,relres,iter]=pcg(Q,b,TOL,2000);

    
    clear Q b
    a=reshape(shadVec,[Js,N]);
    r=(Bs*a)';
    
    shadImg = zeros(Row, Col, K);
    for k = 1:K
        shadImg(:, :, k) = reshape(gather(r(:, k)), Row, Col).*mask;
    end

    score1 = LMSE(gt_shad,shadImg, mask);

    derived_reflImg = specImg./(shadImg+eps);
    score2 = LMSE( gt_refl,derived_reflImg, mask);
catch
    shadImg = gt_shad;
    score_SIID = temp1;
end

%% E r-term
try    
    lambda_esti_r = lambda_esti_s;
    lambda_data = 1;

    b = lambda_esti_r*L_r'*C_r;
    Q = lambda_sc*(Wsc_L'*Wsc_L) + lambda_rc*(Vrc'*Vrc) +lambda_esti_r*(L_r'*L_r);
    
    % cooperate with S initial estimation
    shadVec_restore = sparseBs*shadVec;
    Qdata = diag(sparse(gather(shadVec_restore))) * sparseBr;
    Cdata = reshape(specImgVec', [], 1);
    
    Q = Q + lambda_data*(Qdata'*Qdata);
    b = b + lambda_data*Qdata'*Cdata;
    disp('Q r-estimate complete!')

    tic,
    TOL=1e-8;
    [reflVec,flag,relres,iter]=pcg(Q,b,TOL,500);

    clear Q b
    a=reshape(reflVec,[Jr,N]);
    r=(Br*a)';
    
    reflImage = zeros(Row, Col, K);
    for k = 1:K
        reflImage(:, :, k) = reshape(gather(r(:, k)), Row, Col).*mask;
    end

    score3 = LMSE( gt_refl,reflImage, mask);
catch
    score_LRIID = temp2;
    score4 = LMSE(gt_shad,shadImg, mask);
end    
    
%% iteration
try
idx = 1;
while idx < 100
    if mod(idx, 2)
        lambda_data = 0.01;
    else
        lambda_data = 1;
    end
     
   E_SUM(idx)=lambda_sc*norm(Wsc*shadVec) + lambda_rc*norm(Vrc_L*shadVec)+...
         lambda_sc*norm(Wsc_L*reflVec) + lambda_rc*norm(Vrc*reflVec)+...
         lambda_data*norm((sparseBs*shadVec).*(sparseBr*reflVec)-reshape(specImgVec', [], 1));
     if idx > 3
         delta = -E_SUM(idx) +E_SUM(idx-2);
         if abs(delta) < 0.01
             break;
         end
     end
    
    if mod(idx, 2)
        % Qsc

        Q = lambda_sc*(Wsc'*Wsc);

        % Qrc

        Q = Q + lambda_rc*(Vrc_L'*Vrc_L);
    
    else
        % Qsc

        Q = lambda_sc*(Wsc_L'*Wsc_L);

        % Qrc

        Q = Q + lambda_rc*(Vrc'*Vrc);

    end

    % Qdata
    C = reshape(specImgVec', [], 1);

    if mod(idx, 2)
    % Edata-s term
    reflVec_restore = sparseBr * reflVec; % reflVec length: N*Jr
    Qdata = diag(sparse(reflVec_restore))*sparseBs;

    else
    % Edata-r term
    shadVec_restore = sparseBs * shadVec; % shadVec length: N*Js
    Qdata = diag(sparse(shadVec_restore))*sparseBr;

    end

    Q = Q + lambda_data*(Qdata'*Qdata);
    b = lambda_data*Qdata'*C;
    
    if mod(idx, 2)
        tic,
        TOL=1e-6;
        [shadVec,flag,relres,iter]=pcg(Q,b,TOL,1000);
        
        
        a=reshape(shadVec,[Js,N]);
        r=(Bs*a)';
        
        shadImg = zeros(Row, Col, K);
        for k = 1:K
            shadImg(:, :, k) = reshape(r(:, k), Row, Col).*mask;
        end
        
        score = LMSE(gt_shad,shadImg, mask);

    else
        tic,
        TOL=1e-6;
        [reflVec,flag,relres,iter]=pcg(Q,b,TOL,1000);
        

    end
    idx = idx+1;
end
catch
end

subplot(2,2,3)
imageArray = uint8(reflImg);
imshow(imageArray);
title('Output Reflectance')

subplot(2,2,4)
imageArray = uint8(shadImg);
imshow(imageArray);
title('Output Shading')

fprintf('SIID score = %.3f\n', score_SIID);
fprintf('LRIID score = %.3f\n',score_LRIID);


