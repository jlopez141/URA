set terminal pdf enhanced font ",12"
set output "plot_pdf.pdf"

set xlabel "r"
set ylabel "pdf_{i}"

plot for [t = 2:8] "contr_tot_gdr12.dat" u 1:t w l title "pair ".(t-1)