%%Exemple simple de noeud suscriber matlab avec une structure permettant la
%%genération de code c++ avec codegen

function myNode  

rosshutdown; %suppression du dernier node specifique matlab
rosinit('10.42.0.133'); %création du node spécifique matlab, adresse ip du master
sub = rossubscriber('/point','geometry_msgs/Point',...
    @callback,...
    'DataFormat','struct')

fprintf('Created %s suscriber\n',sub.TopicName);
while(1)
    fprintf('Node is alive..\n');

    pausesub.LatestMessage(3);
end

end %myNode

%Suscriber callback function
function callback(~,msg)
fprintf("(X,Y,Z : (%f,%f,%f)\n",msg.X,msg.Y,msg.Z);
pause(1)
sub.LatestMessage
end