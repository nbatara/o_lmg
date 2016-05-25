
for i=1:21
figure, surf(x*1e9,y*1e9,squeeze(generation_rate(1,1,:,:,1,i))','EdgeColor','none'), view(2), set(gca,'FontSize',fs),colorbar,axis tight
name=['generation_rate_', num2str(round(pabs1.lambda(i)*1e9)) '_nm']
savefig(name)
print(name,'-dpng')
end