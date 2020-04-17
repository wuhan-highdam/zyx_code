clear
fn_av = [];
ft_av = [];
for f=0:50;
    filename1=['frame_all\CNS_frame_',num2str(f),'.dat'];
    inf=dlmread(filename1);
    fn_av = [fn_av;mean(inf(:,1))];
    ft_av = [ft_av;mean(inf(:,5))];
end
