function detect_array=ellipse_interface_scatter(a,b,a_ell,b_ell,quadrant,center_x,center_y,points)
%enter the original x-y coordinate data
x_vec_1=linspace(center_x+a_ell,center_x,points)*2/a;
y_vec_1=( ( 2*b_ell*( 1-( (a*x_vec_1-2*center_x)/(2*a_ell) ).^2 ).^(1/2) )+2*center_y )/b;
detect_array1=[x_vec_1' y_vec_1'];
x_vec_2=linspace(center_x,center_x-a_ell,points)*2/a;
y_vec_2=( ( 2*b_ell*( 1-( (a*x_vec_2-2*center_x)/(2*a_ell) ).^2 ).^(1/2) )+2*center_y )/b;
detect_array2=[x_vec_2' y_vec_2'];
x_vec_3=linspace(center_x-a_ell,center_x,points)*2/a;
y_vec_3=( (-( 2*b_ell*( 1-( (a*x_vec_3-2*center_x)/(2*a_ell) ).^2 ).^(1/2) ) )+2*center_y )/b;
detect_array3=[x_vec_3' y_vec_3'];
x_vec_4=linspace(center_x,center_x+a_ell,points)*2/a;
y_vec_4=( (-2*b_ell*( 1-( (a*x_vec_4-2*center_x)/(2*a_ell) ).^2 ).^(1/2) )+2*center_y )/b;
detect_array4=[x_vec_4' y_vec_4'];
if quadrant==1
    detect_array=detect_array1;
elseif quadrant==2
    detect_array=detect_array2;
elseif quadrant==3
    detect_array=detect_array3;
elseif quadrant==4
    detect_array=detect_array4;
elseif quadrant==12
    detect_array=[detect_array1;detect_array2];
elseif quadrant==123
    detect_array=[detect_array1;detect_array2;detect_array3];
elseif quadrant==1234
    detect_array=[detect_array1;detect_array2;detect_array3;detect_array4];
elseif quadrant==23
    detect_array=[detect_array2;detect_array3];
elseif quadrant==234
    detect_array=[detect_array2;detect_array3;detect_array4];
elseif quadrant==34
    detect_array=[detect_array3;detect_array4];
end
detect_array=real(detect_array);