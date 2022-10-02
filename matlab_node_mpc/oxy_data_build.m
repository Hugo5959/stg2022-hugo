load("donoxori.mat");
load("donox.mat")
donox1 = [];
donox2 = [];
donox3 = [];
donox4 = [];
donox5 = [];
donox6 = [];
donox7 = [];
donox8 = [];
donox9 = [];
donox10 = [];
donox11 = [];
donox12 = [];
donox13 = [];
donox14 = [];
donox15 = [];
donox16 = [];
donox17 = [];
t = 4783;

for i=1:length(donoxori)
    donox1(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end

for i=1:length(donoxori)
    donox2(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end

for i=1:length(donoxori)
    donox3(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end 


for i=1:length(donoxori)
    donox4(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end

for i=1:length(donoxori)
    donox5(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);

end

for i=1:length(donoxori)
    donox6(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end
for i=1:length(donoxori)
    donox7(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t +(1/0.7);
end
for i=1:length(donoxori)
    donox8(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + 0.8;
 
end

for i=1:length(donoxori)
    donox9(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end

for i=1:length(donoxori)
    donox10(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);

end


for i=1:length(donoxori)
    donox11(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end
for i=1:length(donoxori)
    donox12(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end

for i=1:length(donoxori)
    donox13(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end
for i=1:length(donoxori)
    donox14(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
 
end

for i=1:length(donoxori)
    donox15(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end

for i=1:length(donoxori)
    donox16(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);
end

for i=1:length(donoxori)
    donox17(end+1) = donoxori(i) + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;
    t = t + (1/0.7);

end

donox_dyn = [donox;donox1;donox2;donox3;donox4;donox5;donox6;donox7;donox8;donox9;donox10;donox11;donox12;donox13;donox14;donox15;donox16;donox17];






