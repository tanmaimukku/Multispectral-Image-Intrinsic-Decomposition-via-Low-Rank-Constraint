function [resize,alpha,beta,lambda_sc,lambda_rc,lambda_esti_s,lambda_esti_r,lambda_data]=parameters(input)
    
    num = length(input);
    lambda=0;
    if (num<1),	 lambda          = 250;      end
    if (num<2),	 lambda           = 5000;     end
    if (num<3),  lambda            = 3e-4;     end
    if (num<4),  lambda       = 100;      end
    if (num<5),  lambda       = 100;      end
    if (num<6),  lambda   = 0.01;     end
    if (num<7),  lambda   = 0.01;     end
    if (num<8),  lambda     = 1;        end

end
function [] = LRIID_highlight(specImg1,gt_shad1,mask1,varargin1)
    
    [resize,alpha,beta,lambda_sc,lambda_rc,lambda_esti_s,lambda_esti_r,lambda_data] = parameters(varargin);


    gt_refl =  specImg./(gt_shad+eps); 
    [mask, specImg, gt_shad, gt_refl] = RESIZE(resize, mask, specImg, gt_shad, gt_refl);
    [row, col, band] = size(specImg);
    num = row*col;

    specImgVec = reshape(specImg, [], band);
    specImg = specImg/(max(specImg(:))- min(specImg(:)));
    gt_shad = gt_shad/(max(gt_shad(:))- min(gt_shad(:)));
    Jr = 8;
    load 450to700_118chan_basis.mat SPECBASIS
    coeff = pca(SPECBASIS(:, 1:3:end));
    Br = coeff(:, 1:Jr);
    %-----------Bs---------------
    Js = 1;
    load illumination.mat illuminationREDHEAD
    Bs = illuminationREDHEAD(1:3:end)';
    Bs = Bs/norm(Bs);

    %% Weight Map
    [Wr, Wd] = weightMap(Normalized(specImg), alpha, beta);

    %% Weight Matrices
    NeighSize = 2;
    rowIndex = zeros(1,NeighSize*num*band*2); rowIndex(1:2) = [1,1]; 
    colIndex = zeros(1,NeighSize*num*band*2);

    WscVal = zeros(1,NeighSize*num*band*2);
    Wsc_LVal = zeros(1,NeighSize*num*band*2);
    Vrc_LVal = zeros(1,NeighSize*num*band*2);
    VrcVal = zeros(1,NeighSize*num*band*2);

    idx = 1;
    for Neighbor=1:NeighSize
        for i = 1:num
            coor = [i+1,i+row];
            p = i;
            if (coor(Neighbor)<num) && (coor(Neighbor)>0)
                q = coor(Neighbor);
                if Neighbor ==1
                    w = Wd(p);
                else
                    w = Wr(p);
                end             
                v = 1-w;

                rowIndex((idx-1)*2*band+1:2:idx*2*band-1) = (idx-1)*band+1:idx*band;
                rowIndex((idx-1)*2*band+2:2:idx*2*band)   = (idx-1)*band+1:idx*band;
                colIndex((idx-1)*2*band+1:2:idx*2*band-1) = (p-1)*band+1:p*band;
                colIndex((idx-1)*2*band+2:2:idx*2*band)   = (q-1)*band+1:q*band;

                WscVal((idx-1)*2*band+1:2:idx*2*band-1)   = w;
                WscVal((idx-1)*2*band+2:2:idx*2*band)     = -w;
                Wsc_LVal((idx-1)*2*band+1:2:idx*2*band-1) = w*specImgVec(q,:);
                Wsc_LVal((idx-1)*2*band+2:2:idx*2*band)   = -w*specImgVec(p,:);
                Vrc_LVal((idx-1)*2*band+1:2:idx*2*band-1) = v*specImgVec(q,:);
                Vrc_LVal((idx-1)*2*band+2:2:idx*2*band)   = -v*specImgVec(p,:);
                VrcVal((idx-1)*2*band+1:2:idx*2*band-1)   = v;
                VrcVal((idx-1)*2*band+2:2:idx*2*band)     = -v;

                idx = idx+1;
            end
        end
    end

    % Angela: ���ֱ߽��û���������򣬶��ಿ���ÿգ�sparse��
    zeroMark = find(rowIndex == 0);
    rowIndex(zeroMark) = [];    colIndex(zeroMark) = [];
    WscVal(zeroMark) = [];      VrcVal(zeroMark) = [];	
    Wsc_LVal(zeroMark) = [];	Vrc_LVal(zeroMark) = [];

    %---------------------sparseBs-----------------------

    Bs_temp = sparse(Bs);
    Bs_temp = repmat({Bs_temp}, num, 1);
    sparseBs = blkdiag(Bs_temp{:});

    % S estimate
    Wsc = sparse(rowIndex, colIndex, WscVal, num*NeighSize*band, num*band)*sparseBs;
    Vrc_L = sparse(rowIndex, colIndex, Vrc_LVal, num*NeighSize*band, num*band)*sparseBs;
    % Angela: ����Ҫ��ʾԼ��ϡ������С����Ϊ���������ÿ�

    % Generic_s
    %---------------------------
    C = reshape(specImgVec', [], 1);

    %---------------------sparseBr-----------------------
    Br_temp = sparse(Br);
    Br_temp = repmat({Br_temp}, num, 1);
    sparseBr = blkdiag(Br_temp{:});

    % R estimate
    Wsc_L = sparse(rowIndex, colIndex, Wsc_LVal, num*NeighSize*band, num*band)*sparseBr;
    Vrc = sparse(rowIndex, colIndex, VrcVal, num*NeighSize*band, num*band)*sparseBr;

    %% Estimate s

    b = lambda_esti_s*sparseBs'*C;
    Q = lambda_sc*(Wsc'*Wsc) + lambda_rc*(Vrc_L'*Vrc_L) +lambda_esti_s*(sparseBs'*sparseBs);%+ lambda_gc*(Gc'*Gc);
    disp('Q s-estimate complete!')


    % ----------------conjugate gradient descend---------------
    TOL=1e-8;
    [shadVec,~,~,~]=pcg(gpuArray(Q),gpuArray(b),TOL,2000);
    shadVec = gather(shadVec);

    a=reshape(shadVec,[Js,num]);
    r=(Bs*a)';

    shadImg = zeros(row, col, band);
    for k = 1:band
        shadImg(:, :, k) = reshape(r(:, k), row, col).*mask;
    end

    score1 = LMSE(gt_shad,shadImg, mask)

    derived_reflImg = specImg./(shadImg+eps);
    score2 = LMSE( gt_refl,derived_reflImg, mask)

    %% E r-term



    b = lambda_esti_r*sparseBr'*C;
    Q = lambda_sc*(Wsc_L'*Wsc_L) + lambda_rc*(Vrc'*Vrc) +lambda_esti_r*(sparseBr'*sparseBr);

    % cooperate with S initial estimation
    Qdata = diag(sparse(sparseBs*shadVec)) * sparseBr;

    Q = Q + lambda_data*(Qdata'*Qdata);
    b = b + lambda_data*Qdata'*C;
    disp('Q r-estimate complete!')

    TOL=1e-8;
    [reflVec,~,~,~]=pcg(gpuArray(Q),gpuArray(b),TOL,500);
    reflVec = gather(reflVec);

    a=reshape(reflVec,[Jr,num]);
    r=(Br*a)';

    reflImage = zeros(row, col, band);
    for k = 1:band
        reflImage(:, :, k) = reshape(r(:, k), row, col).*mask;
    end

    score3 = LMSE( gt_refl,reflImage, mask)

    %% iteration
    idx = 1;
    while idx < 10
        if mod(idx, 2)
            lambda_data = 0.01;
        else
            lambda_data = 1;
        end

        E_SUM(idx)=lambda_sc*norm(Wsc*shadVec) + lambda_rc*norm(Vrc_L*shadVec)+...
             lambda_sc*norm(Wsc_L*reflVec) + lambda_rc*norm(Vrc*reflVec)+...
             lambda_data*norm((sparseBs*shadVec).*(sparseBr*reflVec)-C); 
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
        if mod(idx, 2)
            % Edata-s term
            Qdata = diag(sparse(sparseBr*reflVec))*sparseBs;
        else
            % Edata-r term
            Qdata = diag(sparse(sparseBs*shadVec))*sparseBr;
        end

        Q = Q + lambda_data*(Qdata'*Qdata);
        b = lambda_data*Qdata'*C;

        if mod(idx, 2)
            TOL=1e-6;
            [shadVec,~,~,~]=pcg(gpuArray(Q),gpuArray(b),TOL,1000);
            shadVec = gather(shadVec);

            a=reshape(shadVec,[Js,num]);
            r=(Bs*a)';

            shadImg = reshape(r,row,col,band).*repmat(mask,1,1,band);

            score = LMSE(gt_shad,shadImg, mask) 

        else
            TOL=1e-6;
            [reflVec,~,~,~]=pcg(gpuArray(Q),gpuArray(b),TOL,1000);
            reflVec = gather(reflVec);
        end
        idx = idx+1;
    end
end

