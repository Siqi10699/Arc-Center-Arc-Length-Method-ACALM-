function Truss_planar_3_members
%% Uncomment one line at a time to reproduce benchmark Ⅰ
E1=0.75;max_step=17;lambda_initial=0.1;max_iteration=5;fname = "input_Truss_2D_3members_model1.txt";
% E1=0.50;max_step=23;lambda_initial=0.1;max_iteration=5;fname = "input_Truss_2D_3members_model2.txt";
% 
% E1=0.75;max_step=17;lambda_initial=0.1;max_iteration=20;fname = "input_Truss_2D_3members_model1.txt";
% E1=0.50;max_step=23;lambda_initial=0.1;max_iteration=20;fname = "input_Truss_2D_3members_model2.txt";
% 
% E1=0.75;max_step=17;lambda_initial=0.1;max_iteration=50;fname = "input_Truss_2D_3members_model1.txt";
% E1=0.50;max_step=23;lambda_initial=0.1;max_iteration=50;fname = "input_Truss_2D_3members_model2.txt";

% E1=0.75;max_step=8;lambda_initial=0.3;max_iteration=5;fname = "input_Truss_2D_3members_model1.txt";
% E1=0.50;max_step=19;lambda_initial=0.3;max_iteration=5;fname = "input_Truss_2D_3members_model2.txt";

% E1=0.75;max_step=11;lambda_initial=0.5;max_iteration=5;fname = "input_Truss_2D_3members_model1.txt";
% E1=0.50;max_step=16;lambda_initial=0.5;max_iteration=5;fname = "input_Truss_2D_3members_model2.txt";
%% Data preparation
[ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r, dof_force, Q, maxloadSteps, loadincr, outputlist] = processfile(fname); %#ok<ASGLU>
% Set uw_or_vw_draw to 1 if 2D projection plots are required
uw_or_vw_draw=0;
u_0=zeros(size(Q));lambda_0=0;clin_sph=1;
u_i=u_0;lambda_i=lambda_0+lambda_initial;i_save=1;non_convergence_count=0;multiply_factor=1.25;reduction_factor=4;
u_iter_save=zeros(size(u_0,1),max_step*max_iteration);error_iter_save=u_iter_save;lambda_iter_save=zeros(max_step*max_iteration,1);
sign_A_save=zeros(max_step,1);per_step_iteration_number_save=sign_A_save;arc_length_save=sign_A_save;
arc_center=zeros(size(u_0(assy4r),1)+1,max_step);lambda_step_save=sign_A_save;fail_step=0;fail_step_save=[];
u_step_save=zeros(size(u_0,1),max_step);error_step_save=u_step_save;detect_array_succ=[];detect_array_fail=detect_array_succ;
sqrt_Q=sqrt(Q'*Q);complex_roots=0;bf=0;K_tangent = zeros(neq,neq);R = zeros(neq,1);tolerance=1e-6;
% Finite element assembly for the initial configuration
if(ndim == 2)
    if(ndof == 2) % Truss element
        for e = 1:nelem
            [Klocal, Flocal] = Truss_2D_model1(elemData, elemConn, e, coords, u_i, bf);
            K_tangent = Assembly_Matrix(K_tangent,Klocal,LM,e);
            R = Assembly_Vector(R,Flocal,LM,e);
        end
    else % Beam element
        for e = 1:nelem
            [Klocal, Flocal] = GeomExactBeam_2D(elemData, elemConn, e, coords, u_i, bf);
            K_tangent = Assembly_Matrix(K_tangent,Klocal,LM,e);
            R = Assembly_Vector(R,Flocal,LM,e);
        end
    end
else
    if(ndof == 3) % Truss element
        for e = 1:nelem
            [Klocal, Flocal] = Truss_3D_model2(elemData, elemConn, e, coords, u_i, bf);
            K_tangent = Assembly_Matrix(K_tangent,Klocal,LM,e);
            R = Assembly_Vector(R,Flocal,LM,e);
        end
    end
end
R=R+lambda_initial*Q;
K_tangent_reduced=K_tangent(assy4r,assy4r);R_reduced=R(assy4r);Q_reduced=Q(assy4r);
%% Iteration procedure
for p=1:max_step
    for i=0:max_iteration
        delta_u_Ri=K_tangent_reduced\R_reduced;
        delta_u_Qi=K_tangent_reduced\Q_reduced;
        % Predictor stage (i=0)
        if i==0
            temp_count=1+non_convergence_count;
            if p==1
                temp_save_for_non_convergence=[u_i(assy4r) delta_u_Qi];
                delta_L=lambda_initial*sqrt(delta_u_Qi'*delta_u_Qi+Q_reduced'*Q_reduced*clin_sph);delta_u_i=lambda_initial*delta_u_Qi;
                u_i_plus1=delta_u_i;lambda_i_plus1=lambda_initial;
            elseif p==temp_count
                delta_u_Qi=temp_save_for_non_convergence(:,2);
                delta_L=delta_L/reduction_factor;lambda_initial=lambda_initial/reduction_factor;delta_u_i=lambda_initial*delta_u_Qi;
                u_i_plus1=delta_u_i;lambda_i_plus1=lambda_initial;
            else
                if last_step_success==1
                    temp_save_for_non_convergence=[u_i(assy4r) delta_u_Qi];lambda_i_save_for_non_convergence=lambda_i;
                    delta_L=delta_L*multiply_factor;
                elseif last_step_success==0
                    u_i(assy4r)=temp_save_for_non_convergence(:,1);lambda_i=lambda_i_save_for_non_convergence;
                    delta_u_Qi=temp_save_for_non_convergence(:,2);
                    delta_L=delta_L/reduction_factor;
                end
                delta_lambda_i=delta_L/sqrt(delta_u_Qi'*delta_u_Qi+Q_reduced'*Q_reduced*clin_sph);
                u_i_plus1_temp=u_i(assy4r)+[delta_lambda_i*delta_u_Qi -delta_lambda_i*delta_u_Qi];
                lambda_i_plus1_temp=[lambda_i+delta_lambda_i lambda_i-delta_lambda_i];
                u_lambda_i_temp=[u_i_plus1_temp;lambda_i_plus1_temp*sqrt_Q];
                arc_center_temp=arc_center(:,p-temp_count);
                % arc center selection criterion
                if norm(arc_center_temp-u_lambda_i_temp(:,1))>norm(arc_center_temp-u_lambda_i_temp(:,2))
                    u_i_plus1=u_i_plus1_temp(:,1);
                    lambda_i_plus1=lambda_i_plus1_temp(1);
                else
                    u_i_plus1=u_i_plus1_temp(:,2);
                    lambda_i_plus1=lambda_i_plus1_temp(2);
                end
            end
            u_p_start=u_i_plus1;lambda_p_start=lambda_i_plus1;
            arc_center(:,p)=[u_i_plus1;sqrt_Q*lambda_i_plus1];
        % Corrector stage (i=1,2,3,...)
        else
            vi=u_i(assy4r)-u_p_start+delta_u_Ri;
            wi=lambda_i-lambda_p_start;
            if clin_sph==0
                a_eq=delta_u_Qi'*delta_u_Qi;
                b_eq=delta_u_Qi'*vi;
                c_eq=vi'*vi-delta_L^2;
            elseif clin_sph==1
                a_eq=delta_u_Qi'*delta_u_Qi+Q_reduced'*Q_reduced;
                b_eq=delta_u_Qi'*vi+wi*(Q_reduced'*Q_reduced);
                c_eq=vi'*vi-delta_L^2+wi^2*(Q_reduced'*Q_reduced);
            end
            if (b_eq^2-a_eq*c_eq)>0
                delta_lambda_i=[(-b_eq-sqrt(b_eq^2-a_eq*c_eq))/a_eq,(-b_eq+sqrt(b_eq^2-a_eq*c_eq))/a_eq];
                delta_u_i=delta_u_Ri+delta_u_Qi*delta_lambda_i;
                u_i_plus1=u_i(assy4r)+delta_u_i;
                lambda_i_plus1=lambda_i+delta_lambda_i;
                % root selection criterion
                if p==1||p==temp_count
                    max_lambda_i_plus1=max(lambda_i_plus1);
                    index=find(lambda_i_plus1==max_lambda_i_plus1);
                    u_i_plus1=u_i_plus1(:,index);lambda_i_plus1=lambda_i_plus1(index);
                else
                    u_step_tempt=u_step_save(assy4r,p-temp_count);lambda_step_tempt=lambda_step_save(p-temp_count);
                    temp_two_roots_distance=[norm([u_i_plus1(:,1)-u_step_tempt;lambda_i_plus1(1)*sqrt_Q-lambda_step_tempt]) norm([u_i_plus1(:,2)-u_step_tempt;lambda_i_plus1(2)*sqrt_Q-lambda_step_tempt])];
                    index=find(temp_two_roots_distance==max(temp_two_roots_distance));
                    u_i_plus1=u_i_plus1(:,index);
                    lambda_i_plus1=lambda_i_plus1(index);
                end
            else
                complex_roots=1;break;
            end
        end

        u_i(assy4r) = u_i_plus1;K_tangent(1:end,1:end) = 0.0;R(1:end) = 0.0;
        u_iter_save(:,i_save)=u_i;
        lambda_iter_save(i_save)=lambda_i_plus1*sqrt_Q;
        % Finite element assembly for the next iteration/step
        if(ndim == 2)
            if(ndof == 2) % Truss element
                for e = 1:nelem
                    [Klocal, Flocal] = Truss_2D_model1(elemData, elemConn, e, coords, u_i, bf);
                    K_tangent = Assembly_Matrix(K_tangent,Klocal,LM,e);
                    R = Assembly_Vector(R,Flocal,LM,e);
                end
            else % Beam element
                for e = 1:nelem
                    [Klocal, Flocal] = GeomExactBeam_2D(elemData, elemConn, e, coords, u_i, bf);
                    K_tangent = Assembly_Matrix(K_tangent,Klocal,LM,e);
                    R = Assembly_Vector(R,Flocal,LM,e);
                end
            end
        else
            if(ndof == 3) % Truss element
                for e = 1:nelem
                    [Klocal, Flocal] = Truss_3D_model2(elemData, elemConn, e, coords, u_i, bf);
                    K_tangent = Assembly_Matrix(K_tangent,Klocal,LM,e);
                    R = Assembly_Vector(R,Flocal,LM,e);
                end
            end
        end
        R = R + lambda_i_plus1*Q;
        K_tangent_reduced=K_tangent(assy4r,assy4r);R_reduced=R(assy4r);Q_reduced=Q(assy4r);
        
        error_iter_save(:,i_save)=R;
        i_save=i_save+1;
        max_R=max(abs(R_reduced));
        fprintf("   Maximum_error/tolerance (%dth iter, %dth step): %10.9f\n",i,p,max_R/tolerance);
        if max_R<tolerance&&i>0
            fprintf("    Final true load:%e      Step %d convergence achieved!!!\n\n",lambda_i_plus1,p);
            lambda_i=lambda_i_plus1;u_i(assy4r)=u_i_plus1;
            break;
        end
        lambda_i=lambda_i_plus1;u_i(assy4r)=u_i_plus1;
    end
    i_last=i;u_step_save(:,p)=u_i;lambda_step_save(p)=lambda_i*sqrt_Q;error_step_save(:,p)=error_iter_save(:,i_save-1);
    arc_length_save(p)=delta_L;per_step_iteration_number_save(p)=i_last+1;
    last_step_success=i_last<=max_iteration&&max_R<tolerance;
    if last_step_success==1
        non_convergence_count=0;
        detect_array_tem=sphere_interface_scatter(u_p_start(2),u_p_start(3),sqrt_Q*lambda_p_start,delta_L,30);
        detect_array_succ=[detect_array_succ;detect_array_tem]; %#ok<AGROW>
    elseif last_step_success==0
        fail_step=fail_step+1;fail_step_save=[fail_step_save;p]; %#ok<AGROW> 
        non_convergence_count=non_convergence_count+1;
        if complex_roots==0
            fprintf("   Step %d convergence failed due to exccessive iterations--\n\n",p);
        else
            fprintf("   Step %d convergence failed due to complex roots--\n\n",p);complex_roots=0;
        end
        detect_array_tem=sphere_interface_scatter(u_p_start(2),u_p_start(3),sqrt_Q*lambda_p_start,delta_L,30);
        detect_array_fail=[detect_array_fail;detect_array_tem]; %#ok<AGROW>
    end
end
fprintf("__________________________________________________________________________________________________________________\n");
fprintf("   fail step:\n");disp(fail_step_save');fprintf("   Total number of fail step: %d\n",fail_step);
u_step_save=u_step_save(:,1:p);lambda_step_save=lambda_step_save(1:p);
keep_columns=true(size(u_step_save,2),1);keep_columns(fail_step_save)=false;
u_step_success_reduced=u_step_save(assy4r,keep_columns);lambda_step_success=lambda_step_save(keep_columns);
fprintf("   Total iteration: %d\n",sum(per_step_iteration_number_save));
fprintf("   Average iter per step (including 0th iter and all failed iter): %3.2f\n",sum(per_step_iteration_number_save)/max_step);
per_step_iteration_number_save(fail_step_save)=0;
fprintf("   Average iter per sucessful step (including 0th iter and all failed iter): %3.2f\n",sum(per_step_iteration_number_save)/(max_step-fail_step));
%% plotting
hold on;axis equal;
scatter3(-detect_array_succ(:,1),-detect_array_succ(:,2),detect_array_succ(:,3),0.001,'blue','.');
if ~isempty(detect_array_fail)
    scatter3(-detect_array_fail(:,1),-detect_array_fail(:,2),detect_array_fail(:,3),0.001,[0.6,0.6,0.6],'.');
end
scatter3(-u_step_success_reduced(2,:)',-u_step_success_reduced(3,:)',lambda_step_success,60,'r','^','filled');
h=sqrt(3)/2;v1 = linspace(0, 2.25*h, 2000);L13 = sqrt(0.25 + (h - v1).^2);
P = -2 * (L13 - 1) .* (h - v1) ./ L13;E1_vals =E1;
for i = 1:length(E1_vals)
    E1 = E1_vals(i);
    v2 = v1 + P / E1;
    scatter3(v1,v2,P,2,'black','o','filled');
end
grid off;xlabel('v1');ylabel('v2');zlabel('Load');
set(groot, 'defaultFigurePosition', [640, 291.5, 747.5, 427]); view(2.273701831961113e+02,17.230624466141293);
xlim([-0.3 2.3]);ylim([-0.1 3.3]);zlim([-0.8 1]);view(90,0);
if uw_or_vw_draw
    figure;xlim([0 3]);ylim([-0.6 0.6]);hold on;set(gca,'color','none');
    scatter(-u_step_success_reduced(2,:)',lambda_step_success,60,'r','^','filled');
    scatter(-u_step_success_reduced(3,:)',lambda_step_success,60,'r','^','filled');
end