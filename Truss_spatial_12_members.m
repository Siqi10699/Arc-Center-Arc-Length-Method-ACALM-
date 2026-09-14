function Truss_spatial_12_members
%% Uncomment one line at a time to reproduce benchmark Ⅱ
multiply_factor=1.25;reduction_factor=4;fig_plot="Fig.9 (b)";max_step=24;
% multiply_factor=1.25;reduction_factor=4;fig_plot="Fig.9 (c)";max_step=24;
% multiply_factor=1.25;reduction_factor=4;fig_plot="Fig.9 (d)";max_step=24;

% multiply_factor=1.25;reduction_factor=3;fig_plot="Fig.10";max_step=21;
% multiply_factor=1.25;reduction_factor=6;fig_plot="Fig.10";max_step=22;
% multiply_factor=1.35;reduction_factor=6;fig_plot="Fig.10";max_step=18;
% multiply_factor=1.15;reduction_factor=3;fig_plot="Fig.10";max_step=21;
%% Data preparation
fname = "input_Truss_3D_12members.txt";
[ndim, ndof, nnode, nelem, coords, elemConn, elemData, LM, neq, assy4r, dof_force, Q, maxloadSteps, loadincr, outputlist] = processfile(fname); %#ok<ASGLU>
lambda_initial=0.1;max_iteration=5;
u_0=zeros(size(Q));lambda_0=0;clin_sph=1;
u_i=u_0;lambda_i=lambda_0+lambda_initial;i_save=1;non_convergence_count=0;
u_iter_save=zeros(size(u_0,1),max_step*max_iteration);error_iter_save=u_iter_save;lambda_iter_save=zeros(max_step*max_iteration,1);
sign_A_save=zeros(max_step,1);per_step_iteration_number_save=sign_A_save;arc_length_save=sign_A_save;
arc_center=zeros(size(u_0(assy4r),1)+1,max_step);lambda_step_save=sign_A_save;fail_step=0;fail_step_save=[];
u_step_save=zeros(size(u_0,1),max_step);error_step_save=u_step_save;
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
    elseif last_step_success==0
        fail_step=fail_step+1;fail_step_save=[fail_step_save;p]; %#ok<AGROW> 
        non_convergence_count=non_convergence_count+1;
        if complex_roots==0
            fprintf("   Step %d convergence failed due to exccessive iterations--\n\n",p);
        else
            fprintf("   Step %d convergence failed due to complex roots--\n\n",p);complex_roots=0;
        end
    end
end
fprintf("__________________________________________________________________________________________________________________\n");
fprintf("   fail step:\n");disp(fail_step_save');fprintf("   Total number of fail step: %d\n",fail_step);
u_step_save=u_step_save(:,1:p);lambda_step_save=lambda_step_save(1:p);
keep_columns=true(size(u_step_save,2),1);keep_columns(fail_step_save)=false;
u_step_success=u_step_save(:,keep_columns);lambda_step_success=lambda_step_save(keep_columns);
u_step_success_reduced__1=u_step_success(10,:);u_step_success_reduced__2=abs(u_step_success(12,:));u_step_success_reduced__3=-u_step_success(15,:);
fprintf("   Total iteration: %d\n",sum(per_step_iteration_number_save));
fprintf("   Average iter per step (including 0th iter and all failed iter): %3.2f\n",sum(per_step_iteration_number_save)/max_step);
per_step_iteration_number_save(fail_step_save)=0;
fprintf("   Average iter per sucessful step (including 0th iter and all failed iter): %3.2f\n",sum(per_step_iteration_number_save)/(max_step-fail_step));
%% plotting
hold on;
set(gca,'color','none');set(groot, 'defaultFigurePosition', [640, 415.5, 541.5, 303]); 
if fig_plot=="Fig.9 (c)"
    scatter(u_step_success_reduced__1, lambda_step_success/sqrt_Q, 60, 'r', '^', 'filled');xlim([0 0.8]); ylim([-0.1 0.1]);
    offset_x = (0.8 - 0) * 0.03;offset_y = (0.1 - (-0.1)) * 0.04;x_text = u_step_success_reduced__1 - offset_x;y_text = (lambda_step_success/sqrt_Q) + offset_y;
    for i = 1:length(u_step_success_reduced__1)
        text(x_text(i), y_text(i), num2str(i),'Color', 'r','FontWeight', 'bold','FontSize', 10,'FontName', 'Times New Roman');
    end
elseif fig_plot=="Fig.9 (b)"
    scatter(u_step_success_reduced__2, lambda_step_success/sqrt_Q, 60, 'r', '^', 'filled');xlim([0 2]); ylim([-0.1 0.1]);
    offset_x = (2 - 0) * 0.03;offset_y = (0.1 - (-0.1)) * 0.04;x_text = u_step_success_reduced__2 - offset_x;y_text = (lambda_step_success/sqrt_Q) + offset_y;
    for i = 1:length(u_step_success_reduced__2)
        text(x_text(i), y_text(i), num2str(i),'Color', 'r','FontWeight', 'bold','FontSize', 10,'FontName', 'Times New Roman');
    end
elseif fig_plot=="Fig.9 (d)"
    scatter(u_step_success_reduced__3([1:7,16:21]), u_step_success_reduced__1([1:7,16:21]), 60, 'r', '^', 'filled');xlim([-0.5 2.5]); ylim([0 0.9]);
    scatter(u_step_success_reduced__3(8:15), u_step_success_reduced__1(8:15), 20, 'r', '^', 'filled');
    offset_x = 3 * 0.03;offset_y = 0.9 * 0.04;x_text = u_step_success_reduced__3 - offset_x;y_text = u_step_success_reduced__1 + offset_y;
    for i = 1:length(u_step_success_reduced__3)
        if ismember(i,8:15)
            text(x_text(i), y_text(i), num2str(i),'Color', 'r','FontWeight', 'bold','FontSize', 6,'FontName', 'Times New Roman');
        else
            text(x_text(i), y_text(i), num2str(i),'Color', 'r','FontWeight', 'bold','FontSize', 10,'FontName', 'Times New Roman');
        end
    end
elseif fig_plot=="Fig.10"
    scatter(u_step_success_reduced__3, u_step_success_reduced__1, 50, 'r', '^', 'filled');xlim([-0.5 2.5]); ylim([0 0.9]);
end