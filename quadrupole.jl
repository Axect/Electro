using PyPlot

struct Params
  q :: Float64 # charge
  q_circ :: Float64 # charge for charged conductor
  r :: Float64 # radius of circle
end

mutable struct Position
  x :: Float64
  y :: Float64
end

function native_quadrupole(q::Float64, P1::Position, P2::Position, P3::Position, P4::Position)
  (x1, y1) = (P1.x, P1.y);
  (x2, y2) = (P2.x, P2.y);
  (x3, y3) = (P3.x, P3.y);
  (x4, y4) = (P4.x, P4.y);
  return (x,y) -> q * (1/sqrt((x-x1)^2 + (y-y1)^2) - 1/sqrt((x-x2)^2 + (y-y2)^2) + 1/sqrt((x-x3)^2 + (y-y3)^2) - 1/sqrt((x-x4)^2 + (y-y4)^2))
end

function image_quadrupole(p::Params, P1::Position, P2::Position, P3::Position, P4::Position, charged=false)
  # x-axis quadrupole
  (q, r) = (p.q, p.r);
  (x1, y1) = (P1.x, P1.y);
  (x2, y2) = (P2.x, P2.y);
  (x3, y3) = (P3.x, P3.y);
  (x4, y4) = (P4.x, P4.y);
  l1 = sqrt(x1^2 + y1^2); # distance of q1 (Center of Circle = (0,0))
  l2 = sqrt(x2^2 + y2^2); # distance of q2
  l3 = sqrt(x3^2 + y3^2); # distance of q3
  l4 = sqrt(x4^2 + y4^2); # distance of q4
  (a1, b1) = r^2/l1 .* (x1/l1, y1/l1); # distance of q1' 
  (a2, b2) = r^2/l2 .* (x2/l2, y2/l2); # distance of q2'
  (a3, b3) = r^2/l3 .* (x3/l3, y3/l3); # distance of q3'
  (a4, b4) = r^2/l4 .* (x4/l4, y4/l4); # distance of q4'
  q01 = q * r/l1; # image charge of q1 (=q1')
  q02 = q * r/l2; # image charge of q2 (=q2')
  q03 = q * r/l3; # image charge of q3 (=q3')
  q04 = q * r/l4; # image charge of q4 (=q4')

  if !charged
    q_circ = 0;
  else
    q_circ = p.q_circ;
  end

  return (x,y) -> q * (1/sqrt((x-x1)^2 + (y-y1)^2) - 1/sqrt((x-x2)^2 + (y-y2)^2)) + 1/sqrt((x-x3)^2 + (y-y3)^2) - 1/sqrt((x-x4)^2 + (y-y4)^2) - q01 * (1/sqrt((x-a1)^2 + (y-b1)^2)) + q02 * (1/sqrt((x-a2)^2 + (y-b2)^2))   - q03 * (1/sqrt((x-a3)^2 + (y-b3)^2)) + q04 * (1/sqrt((x-a4)^2 + (y-b4)^2))+ q_circ * (1/sqrt(x^2 + y^2))
end

function main()
  p1 = Position(3,4);
  p2 = Position(-3,4);
  p3 = Position(-3,-4);
  p4 = Position(3,-4);
  q = 10.0;
  q_circ = 6.0;
  r = 2.0;

  p = Params(q,q_circ,r);

  pot = native_quadrupole(q,p1,p2,p3,p4);
  pot2 = image_quadrupole(p,p1,p2,p3,p4,true);

  x = -10:0.1:10;
  y = -10:0.1:10;

  X = repmat(x', length(y), 1);
  Y = repmat(y, 1, length(x));
  Z = map(pot,X,Y);
  #Z = map(pot2,X,Y);
  
  surf(x,y,Z)
  contour(x,y,Z)
  savefig("native_quadrupole.png", dpi=600)
  #savefig("image_quadrupole.png")
  #savefig("charged_quadrupole.png", dpi=300)
end

main()
