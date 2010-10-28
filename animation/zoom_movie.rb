# zoom_movie.rb: repeatedly runs the fractal program, changing the
#   size, tolerance and filename to produce a "zoom-in" animation.
#
nframes = 1000
(nframes).times do |i|
  ex = 0.2 + (-3.2) * i.to_f / (nframes-1)

  lx = 2.5 * 10**(ex)
  tol = 100/lx

  cmd = "mpirun -np 5 ./fractal_mpi.x fractal/Lx=#{lx} image/tolerance=#{tol} image/file=frame_#{i}.ppm"
  system cmd
  cmd = "convert frame_#{i}.ppm frames/frame_#{i}.png"
  system cmd
  cmd = "/bin/rm frame_#{i}.ppm"
  system cmd
end
