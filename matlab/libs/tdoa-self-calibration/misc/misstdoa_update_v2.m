function solut = misstdoa_update_v2(sol,dz)

m1 = length(sol.rows);
n1 = length(sol.cols);

solut = sol;

Ju1 = 1:m1;
Ju2 = (1:m1)+m1;
Ju3 = (1:m1)+2*m1;
Jv1 = (1:n1) + 3*m1;
Jv2 = (1:n1)+3*m1+1*n1;
Jv3 = (1:n1)+3*m1+2*n1;
Ja  = (1:m1)+3*m1+3*n1;
Jb  = (1:n1)+4*m1+3*n1;
Jo  = (1:n1)+4*m1+4*n1;

solut.u(:,1) = solut.u(:,1)+dz(Ju1);
solut.u(:,2) = solut.u(:,2)+dz(Ju2);
solut.u(:,3) = solut.u(:,3)+dz(Ju3);
solut.v(1,:) = solut.v(1,:)+dz(Jv1)';
solut.v(2,:) = solut.v(2,:)+dz(Jv2)';
solut.v(3,:) = solut.v(3,:)+dz(Jv3)';
solut.a = solut.a + dz(Ja);
solut.b = solut.b + dz(Jb)';
solut.o = solut.o + dz(Jo)';

