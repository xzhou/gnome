function seq = alleleTranslate(seq,oldAlleleLable, newAlleleLable)

old1 = strcat(oldAlleleLable{1,1:end});
new1 = strcat(newAlleleLable{1,1:end});
old2 = strcat(oldAlleleLable{2,1:end});
new2 = strcat(newAlleleLable{2,1:end});

indx1 = (seq == old1);
indx2 = ~indx1;
seq(indx1) = new1(indx1);
seq(indx2) = new2(indx2);
