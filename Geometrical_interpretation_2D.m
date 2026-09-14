function Geometrical_interpretation_2D
%% Data preparation
u_0=0;lambda_0=0;delta_lambda_bar=0.1;max_step=32;max_iteration=5;
clin_sph=1;tolerance=1e-6;u_i=u_0;lambda_i=lambda_0+delta_lambda_bar;i_save=1;detect_array_succ=[];multiply_factor=1.25;reduction_factor=4;
u_iter_save=zeros(size(u_0,1),max_step*max_iteration);error_iter_save=u_iter_save;lambda_iter_save=zeros(max_step*max_iteration,1);
sign_A_save=zeros(max_step,1);per_step_iteration_number_save=sign_A_save;arc_length_save=sign_A_save;detect_array_fail=detect_array_succ;
arc_center=zeros(2,max_step);lambda_step_save=sign_A_save;u_step_save=u_iter_save;error_step_save=u_iter_save;
fail_step=0;fail_step_save=[];Q=5;complex_roots=0;approxi_value=u_i^3-57/8*u_i^2+51/4*u_i;arc_center_succ=[];arc_center_fail=[];
R=approxi_value-Q*lambda_i;sqrt_Q=sqrt(Q'*Q);non_convergence_count=0;lambda_initial=delta_lambda_bar;
%% Iteration procedure
for p=1:max_step
    for i=0:max_iteration
        K_tangent=3*u_i^2-57/4*u_i+51/4;
        delta_u_Ri=-K_tangent\R;
        delta_u_Qi=K_tangent\Q;
        % Predictor stage (i=0)
        if i==0
            temp_count=1+non_convergence_count;
            if p==1
                temp_save_for_non_convergence=[u_i delta_u_Qi];
                delta_L=delta_lambda_bar*sqrt(delta_u_Qi'*delta_u_Qi+Q'*Q*clin_sph);delta_u_i=delta_lambda_bar*delta_u_Qi;
                u_i_plus1=u_0+delta_u_i;lambda_i_plus1=lambda_initial;
            elseif p==temp_count
                delta_u_Qi=temp_save_for_non_convergence(:,2);
                delta_L=delta_L/reduction_factor;lambda_initial=lambda_initial/reduction_factor;delta_u_i=lambda_initial*delta_u_Qi;
                u_i_plus1=delta_u_i;lambda_i_plus1=lambda_initial;
            else
                if last_step_success==1
                    temp_save_for_non_convergence=[u_i delta_u_Qi];lambda_i_save_for_non_convergence=lambda_i;
                    delta_L=delta_L*multiply_factor;
                elseif last_step_success==0
                    u_i=temp_save_for_non_convergence(:,1);lambda_i=lambda_i_save_for_non_convergence;
                    delta_u_Qi=temp_save_for_non_convergence(:,2);
                    delta_L=delta_L/reduction_factor;
                end
                delta_lambda_i=delta_L/sqrt(delta_u_Qi'*delta_u_Qi+Q'*Q*clin_sph);
                u_i_plus1_temp=u_i+[delta_lambda_i*delta_u_Qi -delta_lambda_i*delta_u_Qi];
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
            vi=u_i-u_p_start+delta_u_Ri;
            wi=lambda_i-lambda_p_start;
            if clin_sph==0
                a_eq=delta_u_Qi'*delta_u_Qi;
                b_eq=delta_u_Qi'*vi;
                c_eq=vi'*vi-delta_L^2;
            elseif clin_sph==1
                a_eq=delta_u_Qi'*delta_u_Qi+Q'*Q;
                b_eq=delta_u_Qi'*vi+wi*(Q'*Q);
                c_eq=vi'*vi-delta_L^2+wi^2*(Q'*Q);
            end
            if (b_eq^2-a_eq*c_eq)>0
                delta_lambda_i=[(-b_eq-sqrt(b_eq^2-a_eq*c_eq))/a_eq,(-b_eq+sqrt(b_eq^2-a_eq*c_eq))/a_eq];
                delta_u_i=delta_u_Ri+delta_u_Qi*delta_lambda_i;
                u_i_plus1=u_i+delta_u_i;
                lambda_i_plus1=lambda_i+delta_lambda_i;
                % root selection criterion
                if p==1||p==temp_count
                    max_lambda_i_plus1=max(lambda_i_plus1);
                    index=find(lambda_i_plus1==max_lambda_i_plus1);
                    u_i_plus1=u_i_plus1(:,index);lambda_i_plus1=lambda_i_plus1(index);
                else
                    u_step_tempt=u_step_save(:,p-temp_count);lambda_step_tempt=lambda_step_save(p-temp_count);
                    temp_two_roots_distance=[norm([u_i_plus1(:,1)-u_step_tempt;lambda_i_plus1(1)*sqrt_Q-lambda_step_tempt]) norm([u_i_plus1(:,2)-u_step_tempt;lambda_i_plus1(2)*sqrt_Q-lambda_step_tempt])];
                    index=find(temp_two_roots_distance==max(temp_two_roots_distance));
                    u_i_plus1=u_i_plus1(:,index);
                    lambda_i_plus1=lambda_i_plus1(index);
                end
            else
                complex_roots=1;break;
            end
        end

        u_i = u_i_plus1;u_iter_save(:,i_save)=u_i;
        lambda_iter_save(i_save)=lambda_i_plus1*sqrt_Q;
        approxi_value=u_i_plus1^3-57/8*u_i_plus1^2+51/4*u_i_plus1;
        R=approxi_value-Q*lambda_i_plus1;

        error_iter_save(:,i_save)=R;
        i_save=i_save+1;
        max_R=max(abs(R));
        fprintf("   Maximum_error/tolerance (%dth iter, %dth step): %10.9f\n",i,p,max_R/tolerance);
        if max_R<tolerance&&i>0
            fprintf("    Final true load:%e      Step %d convergence achieved!!!\n\n",lambda_i_plus1,p);
            lambda_i=lambda_i_plus1;u_i=u_i_plus1;
            break;
        end
        lambda_i=lambda_i_plus1;u_i=u_i_plus1;
    end
    i_last=i;u_step_save(:,p)=u_i;lambda_step_save(p)=lambda_i*sqrt_Q;error_step_save(:,p)=error_iter_save(:,i_save-1);
    arc_length_save(p)=delta_L;per_step_iteration_number_save(p)=i_last+1;
    last_step_success=i_last<=max_iteration&&max_R<tolerance;
    if last_step_success==1
        non_convergence_count=0;
        detect_array_tem=ellipse_interface_scatter(2,2,delta_L,delta_L,1234,arc_center(1,p),arc_center(2,p),1500);
        detect_array_succ=[detect_array_succ;detect_array_tem]; %#ok<AGROW>
        arc_center_succ_tem=arc_center(:,p);arc_center_succ=[arc_center_succ_tem arc_center_succ;]; %#ok<AGROW>
    elseif last_step_success==0
        fail_step=fail_step+1;fail_step_save=[fail_step_save;p]; %#ok<AGROW> 
        non_convergence_count=non_convergence_count+1;
        detect_array_tem=ellipse_interface_scatter(2,2,delta_L,delta_L,1234,arc_center(1,p),arc_center(2,p),1500);
        detect_array_fail=[detect_array_fail;detect_array_tem]; %#ok<AGROW>
        arc_center_fail_tem=arc_center(:,p);arc_center_fail=[arc_center_fail_tem arc_center_fail;]; %#ok<AGROW>
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
fprintf("   Total iteration: %d\n",sum(per_step_iteration_number_save));
fprintf("   Average iter per step (including 0th iter and all failed iter): %3.2f\n",sum(per_step_iteration_number_save)/max_step);
%% plotting
hold on;axis equal;set(gca,'color','none');
error_iter_save_final=error_iter_save(:,1:i_save-1); %#ok<NASGU> 
u_step_save=u_step_save(:,1:p);lambda_step_save=lambda_step_save(1:p);
keep_columns=true(size(u_step_save,2),1);keep_columns(fail_step_save)=false;
u_step_save_success=u_step_save(:,keep_columns);lambda_step_save_success=lambda_step_save(keep_columns);
xlabel('X');ylabel('Y');zlabel('Z');
hold on;axis equal;
scatter(detect_array_succ(:,1),detect_array_succ(:,2),0.1,[0,112,192]/256,'.');
if ~isempty(detect_array_fail)
    scatter(detect_array_fail(:,1),detect_array_fail(:,2),0.1,[0.7,0.7,0.7],'.');
end
load('tem_exa.mat'); %#ok<LOAD> 
hold on;plot(tem_exa(:,1),tem_exa(:,2),"k","LineWidth",1,'MarkerSize',1,'Marker','.'); %#ok<USENS> 
scatter(arc_center_succ(1,:),arc_center_succ(2,:),40,[0,112,192]/256,'x');
scatter(arc_center_fail(1,:),arc_center_fail(2,:),40,[0.7,0.7,0.7],'x');
% scatter(u_iter_save_final,lambda_iter_save_final,30,'cyan','^');
scatter((u_step_save_success(1,:))',lambda_step_save_success,50,'r','^','filled');