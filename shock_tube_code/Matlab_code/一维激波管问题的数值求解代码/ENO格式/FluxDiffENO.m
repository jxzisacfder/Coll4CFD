function dF = FluxDiffENO(F,dx,x,direction)
%程序使用统一的正通量迎风格式，负通量时，直接转换成正通量计算
%这里的程序会和边界处理相关
global N;

%根据迎风方向替换实际计算的值
if (direction==1)
    F = F(:,1:end-1);
else
    F = F(:,end:-1:2);
end

%计算三个子模板的差分值
F1 = F(:,1:N+1);
F2 = F(:,2:N+2);
F3 = F(:,3:N+3);
F4 = F(:,4:N+4);
F5 = F(:,5:N+5);

h0 =  1/3*F1 - 7/6*F2 + 11/6*F3;
h1 = -1/6*F2 + 5/6*F3 +  1/3*F4;
h2 =  1/3*F3 + 5/6*F4 -  1/6*F5;

%计算光滑检测子的值
IS0 = 13/12*(F1-2*F2+F3).^2 + 1/4*(F1-4*F2+3*F3).^2;
IS1 = 13/12*(F2-2*F3+F4).^2 + 1/4*(F2-F4).^2;
IS2 = 13/12*(F3-2*F4+F5).^2 + 1/4*(3*F3-4*F4+F5).^2;

%判断使用那个模板确定每个半点处的值
%为了可读性，放弃向量化
%这里有一个小细节，对Euler方程的每一个分量分别选择使用哪个通量并不影响总体的精度和稳定性
for i = 1:N+1
    for n = 1:3
        if (IS0(n,i)<=IS1(n,i) && IS0(n,i)<=IS2(n,i))
            h(n,i) = h0(n,i);
        elseif (IS1(n,i)<=IS0(n,i) && IS1(n,i)<=IS2(n,i))
            h(n,i) = h1(n,i);
        elseif (IS2(n,i)<=IS0(n,i) && IS2(n,i)<=IS1(n,i))
            h(n,i) = h2(n,i);
        end
    end 
end
       
%计算差分   
if (direction==1)
    dF = (h(:,2:end) - h(:,1:end-1))/dx;       
else
    dF = -(h(:,end:-1:2) - h(:,end-1:-1:1))/dx;
end

end



