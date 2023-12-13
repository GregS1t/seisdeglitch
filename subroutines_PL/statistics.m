%%%%%% statistics
%BINNING THE DATA INTO 1 DB SEGMENTS
%Define bin edges
%binedges=linspace(-200,-60,140+1); %(deine the low power and high power, then the difference between the two)
clear pdf_matrix binedges round_PSDs quant_tc whichbintot;
%binedges=linspace(-220,-120,floor(100/6)+1); %(deine the low power and high power, then the difference between the two)
if datat == 'A'
binedges=linspace(-220,-120,floor(100/2)+1); %(deine the low power and high power, then the difference between the two)
else
binedges=linspace(1,100,floor(100/2)+1); %(deine the low power and high power, then the difference between the two)
end
% round each of the power values to a whole number
tot_h=length(P(1,:));
tot_f=length(P(:,1));
for h= 1:tot_h
        round_PSDs(1:tot_f,h)=round(20*log10(sqrt(P(1:tot_f,h))));
end
[numcentral duration]=size(round_PSDs(1:tot_f,1:tot_h));
%determine the number of values in each hour long segment which fall into 
%each 1Db bin. Simultaneously determines which bin each value fell into. 
 for n=1:numcentral
[Np,whichbin]=histc(round_PSDs(n,:),binedges);
%make a matrix out of the quantities pretaining to the number of %estimates which fall into each bin
quant_tc(:,n)=Np;
whichbintot(:,n)=whichbin;
%create the matrix which will be a function of period and power with 
%weight based on the power density function. 
 pdf_matrix(:,n)=quant_tc(:,n)./tot_h;
end

[FF1,BINBIN]=meshgrid(F1,binedges);
stat=pcolor(FF1,BINBIN,pdf_matrix)

set(gca,'Xscale','log')
shading interp
shading flat
if datat == 'A'
ylim(gca,[-220 -120])
else
ylim(gca,[1 100])
end
xlim([5e-4 1.])
ticks=(0:0.1:1);
C=[101,200,182;
104,201,179;
106,202,177;
108,203,174;
111,204,172;
113,204,170;
116,205,167;
118,206,165;
121,207,163;
123,208,160;
125,209,158;
128,210,155;
130,210,153;
133,211,151;
135,212,148;
138,213,146;
140,214,144;
142,215,141;
145,216,139;
147,216,136;
150,217,134;
152,218,132;
155,219,129;
157,220,127;
159,221,125;
162,222,122;
164,223,120;
167,223,117;
169,224,115;
172,225,113;
174,226,110;
176,227,108;
179,228,106;
181,229,103;
184,229,101;
186,230,98;
189,231,96;
191,232,94;
193,233,91;
196,234,89;
198,235,87;
201,235,84;
203,236,82;
205,237,79;
208,238,77;
210,239,75;
213,240,72;
215,241,70;
218,241,68;
220,242,65;
222,243,63;
225,244,60;
227,245,58;
230,246,56;
232,247,53;
235,247,51;
237,248,49;
239,249,46;
242,250,44;
244,251,41;
247,252,39;
249,253,37;
252,253,34;
254,254,32;
255,254,30;
255,252,30;
255,250,29;
255,248,29;
255,246,28;
255,244,28;
255,242,27;
255,240,27;
255,238,26;
255,236,26;
255,234,25;
255,232,25;
255,230,24;
255,228,24;
255,226,23;
255,224,23;
255,222,22;
255,219,22;
255,217,22;
255,215,21;
255,213,21;
255,211,20;
255,209,20;
255,207,19;
255,205,19;
255,203,18;
255,201,18;
255,199,17;
255,197,17;
255,195,16;
255,193,16;
255,191,15;
255,189,15;
255,186,14;
255,184,14;
255,182,13;
255,180,13;
255,178,13;
255,176,12;
255,174,12;
255,172,11;
255,170,11;
255,168,10;
255,166,10;
255,164,9;
255,162,9;
255,160,8;
255,158,8;
255,156,7;
255,154,7;
255,151,6;
255,149,6;
255,147,5;
255,145,5;
255,143,4;
255,141,4;
255,139,4;
255,137,3;
255,135,3;
255,133,2;
255,131,2;
255,129,1;
255,127,1;
255,125,0;
255,123,0;
253,121,0;
251,119,0;
248,117,0;
246,115,0;
244,113,0;
242,111,0;
240,109,0;
238,107,0;
236,105,0;
234,103,0;
232,101,0;
230,99,0;
228,97,0;
226,96,0;
223,94,0;
221,92,0;
219,90,0;
217,88,0;
215,86,0;
213,84,0;
211,82,0;
209,80,0;
207,78,0;
205,76,0;
203,74,0;
201,72,0;
198,70,0;
196,68,0;
194,66,0;
192,64,0;
190,63,0;
188,61,0;
186,59,0;
184,57,0;
182,55,0;
180,53,0;
178,51,0;
176,49,0;
173,47,0;
171,45,0;
169,43,0;
167,41,0;
165,39,0;
163,37,0;
161,35,0;
159,33,0;
157,32,0;
155,30,0;
153,28,0;
150,26,0;
148,24,0;
146,22,0;
144,20,0;
142,18,0;
140,16,0;
138,14,0;
136,12,0;
134,10,0;
132,8,0;
130,6,0;
128,4,0;
125,2,0;
123,0,0;
121,0,0;
119,0,0;
118,0,0;
116,0,0;
114,0,0;
112,0,0;
110,0,0;
108,0,0;
106,0,0;
104,0,0;
102,0,0;
100,0,0;
98,0,0;
96,0,0;
94,0,0;
93,0,0;
91,0,0;
89,0,0;
87,0,0;
85,0,0;
83,0,0;
81,0,0;
79,0,0;
77,0,0;
75,0,0;
73,0,0;
71,0,0;
69,0,0;
67,0,0;
66,0,0;
64,0,0;
62,0,0;
60,0,0;
58,0,0;
56,0,0;
54,0,0;
52,0,0;
50,0,0;
48,0,0;
46,0,0;
44,0,0;
42,0,0;
40,0,0;
39,0,0;
37,0,0;
35,0,0;
33,0,0;
31,0,0;
29,0,0;
27,0,0;
25,0,0;
23,0,0;
21,0,0;
19,0,0;
17,0,0;
15,0,0;
13,0,0;
12,0,0;
10,0,0;
8,0,0;
6,0,0;
4,0,0;
2,0,0;
0,0,0];

set(gca,'XScale','log')
colormap(gca,C/255)
cc=colorbar(gca,'YTick',ticks,'YTickLabel',ticks);
cc.Label.String = 'Probability';
cc.Label.FontSize = 16;
if datat == 'A'
ylabel(gca,'Power [20log_{10} (m/s^2/Hz^{1/2}) ] dB')
else
ylabel(gca,'Power [20log_{10} (DU/Hz^{1/2}) ] dB')
end
xlabel(gca,'Frequency (Hz)')
