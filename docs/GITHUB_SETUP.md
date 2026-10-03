# GitHub setup (new repository)

1. On GitHub: **New repository** → name `dpsc-wls-validation` → Public or Private →
   **do not** add a README, .gitignore or license (they are already here) → Create.

2. In a terminal:

```bash
cd ~/dpsc-wls-validation
git init -b main
git add .
git commit -m "Initial commit: WLS-in-DPSC cross-dataset validation"
git remote add origin git@github.com:mtariqi/dpsc-wls-validation.git
git push -u origin main
```

(Use `https://github.com/mtariqi/dpsc-wls-validation.git` as the remote if you push over HTTPS.)

## Optional: one tracking issue per dataset (GitHub CLI)

```bash
gh label create validation --color 1F5C8B --description "WLS validation" || true
for d in GSE185222 GSE202476 GSE227731
  gh issue create --label validation --title "WLS in DPSCs: $d" --body "See docs/VALIDATION_PLAN.md"
end
```

(The loop above is fish syntax; in bash use `for d in ...; do ...; done`.)
