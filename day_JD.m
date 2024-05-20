function JD = day_JD(year,month,day,hour,minute,second)
% �ó���ʵ�ֹ�����������JD��ת��
% ���÷�Χ1900-2100��

% Check inputs
if nargin < 6
  second = 0;
  if nargin < 5
    minute = 0;
    if nargin < 4
      hour = 0;
      if nargin < 3
        error(message('MATLAB:day_JD:NotEnoughInputs'));
      end  
    end
  end
end
JD = 367*year-floor(1.75*(year+floor((month+9)/12)))+floor(275*month/9)+day+1721013.5+hour/24+minute/(24*60)+second/(24*3600);