int main() {
    printf("backend = %%MPLIB%%\n");
    const char *names = "ESBPNRMO";
    const char *labels[] = {"eps", "sfmin", "base", "precision", "t", "rmin", "rmax"};
    const char codes[] = {'E', 'S', 'B', 'P', 'N', 'U', 'O'};
    for (int i = 0; i < 7; i++) {
        printf("%s = ", labels[i]);
        char c[2];
        c[0] = codes[i];
        c[1] = '\0';
        printnum(Rlamch(c));
        printf("\n");
    }
    return 0;
}
