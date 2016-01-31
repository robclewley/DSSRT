function result = GetColStyleSet(numVars,init_cols,styles,colStep);
% for use by DSSRT
% R. Clewley 2004
res_cols   = zeros(numVars,3);
res_styles = zeros(numVars,1);
rComponent = init_cols(1);
gComponent = init_cols(2);
bComponent = init_cols(3);
styleIx    = 1;
step       = 0;
styleSwitchCount = length(styles);

for varIx=1:numVars
    % switch line styles every styleSwitchCount colours
    step = step + 1;
    if step == styleSwitchCount
        styleIx = mod(styleIx+1,styleSwitchCount);
        if styleIx == 0
            styleIx = styleSwitchCount;
        end
    end
    res_cols(varIx,:)=[rComponent, gComponent, bComponent];
    res_styles(varIx)=styleIx;
    if round((step+1) / 2) == (step+1) / 2 % multiples of 3
        gComponent = gComponent + colStep*0.8; % interdependency 1!
        rComponent = mod(gComponent + colStep*0.7,1);
        if gComponent > 0.9
            gComponent = gComponent-0.9;
        end
        bComponent = mod(bComponent + colStep*0.7,1);
    else
        gComponent = gComponent + colStep;
        if gComponent > 0.8
            gComponent = gComponent-0.8;
            rComponent = rComponent + colStep*0.6;
            if rComponent > 0.85
                rComponent = rComponent-0.85;
                bComponent = bComponent + colStep*0.7;
                if bComponent > 0.8
                    bComponent = bComponent-0.8;
                end
                if abs(bComponent-gComponent)<0.3
                    bComponent = mod(gComponent + 0.3,1);
                end
            end
        end
    end
    diffEnough = false;
    while ~diffEnough & varIx > 1
        doneSearch = true;
        for i = 1:varIx-1
            totdiff = abs(res_cols(i,:)-[rComponent gComponent bComponent]);
            % only worry about close colours if styles are the same
            if sum(totdiff) < 0.5 & res_styles(i) == res_styles(varIx)
                [mindiff mindiffpos] = min(totdiff);
                switch mindiffpos
                    case 1
                        rComponent = mod(rComponent + (2*rand(1)-1)*0.5,1);
                    case 2
                        gComponent = mod(gComponent + (2*rand(1)-1)*0.5,1);
                    case 3
                        bComponent = mod(bComponent + (2*rand(1)-1)*0.5,1);
                end
                doneSearch = false;
                break % for loop
            end
        end
        if doneSearch
            diffEnough = true;
        end
    end
end
result = {res_cols, res_styles};
return