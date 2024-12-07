% load acceleration data (variable named axyz) into the student workspace
load('AccData') % axyz has 34545 rows (time samples) and 3 columns (acceleration along x, y and z axes)

% make sure the function BufferAccData is called for individual time samples, 
% simulating the simulink model running on real time in the student smartphone

N = size(axyz,1); % Time samples
BodyMass = 70; % (kg) mass of the subject from which acceleration data has been saved in AccData.mat
IAAtot = zeros(N,1); % Allocating memory for processed acceleration data
AccBufferComplete = IAAtot; % for AccBufferComplete
EE = IAAtot; % and for Energy Expenditure estimation

% Creating a variable with logical zeros to simulate the event of a user button press
EstimateEEButtonPressed = logical(zeros(size(axyz,1),1));
EstimateEEButtonPressed(50*60:50*60*10) = logical(1); % Simulating the situation whereby the user pressed the button from 1st to 10th minute

for i = 1:N % lover over the 34545 time samples
    [IAAtot(i),AccBufferComplete(i)] = BufferAccData(axyz(i,:)); % Calling the BufferAccData (created in Part 1)
    EE(i) = FromAccToEE(IAAtot(i),AccBufferComplete(i),BodyMass,EstimateEEButtonPressed(i)); % Calling the requested function for the i-th time sample
end

% Creating plot to help students check their solution for the third output argument
ax = axes;
line([1:numel(EE)]/(50*60),EE,'color',[1 1 1]*.8,'linewidth',3,'parent',ax(1))
line([1:numel(IAAtot)]/(50*60),IAAtot,'color','k','parent',ax(1))
xlabel('Time (min)'),ylabel('EE (kcal) and IAA_{tot} (ms^{-1})'), title('Checking: EE should increase with a staircase profile from the 1^{st} to the 10^{th} minute')
legend('EE','IAA_{tot}')