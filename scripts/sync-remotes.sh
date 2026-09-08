#!/usr/bin/env bash
# Bidirectional sync between origin (INSA GitLab) and github (GitHub).
# - Branches only on one side are pushed to the other.
# - When both sides have the branch, only fast-forward pushes are made.
# - Diverged branches are reported without any force-push.

set -euo pipefail

RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'; NC='\033[0m'

echo "Fetching from all remotes..."
git fetch --all --tags --prune --quiet

all_branches=$(git branch -r \
    | grep -E 'remotes/(origin|github)/' \
    | grep -v HEAD \
    | sed 's|.*remotes/\(origin\|github\)/||' \
    | sort -u)

synced=0; skipped=0; conflicts=0

for branch in $all_branches; do
    on_origin=$(git branch -r | grep -c " origin/$branch$" || true)
    on_github=$(git branch -r | grep -c " github/$branch$"  || true)

    if [ "$on_origin" -eq 1 ] && [ "$on_github" -eq 0 ]; then
        echo -e "${GREEN}→ github${NC}: $branch  (origin only)"
        git push github "refs/remotes/origin/$branch:refs/heads/$branch" --quiet
        synced=$((synced + 1))

    elif [ "$on_origin" -eq 0 ] && [ "$on_github" -eq 1 ]; then
        echo -e "${GREEN}→ origin${NC}: $branch  (github only)"
        git push origin "refs/remotes/github/$branch:refs/heads/$branch" --quiet
        synced=$((synced + 1))

    else
        ahead_gh=$(git log --oneline "origin/$branch..github/$branch" 2>/dev/null | wc -l)
        ahead_or=$(git log --oneline "github/$branch..origin/$branch" 2>/dev/null | wc -l)

        if [ "$ahead_gh" -gt 0 ] && [ "$ahead_or" -eq 0 ]; then
            echo -e "${GREEN}→ origin${NC}: $branch  (github ahead by $ahead_gh commit(s))"
            git push origin "refs/remotes/github/$branch:refs/heads/$branch" --quiet
            synced=$((synced + 1))

        elif [ "$ahead_or" -gt 0 ] && [ "$ahead_gh" -eq 0 ]; then
            echo -e "${GREEN}→ github${NC}: $branch  (origin ahead by $ahead_or commit(s))"
            git push github "refs/remotes/origin/$branch:refs/heads/$branch" --quiet
            synced=$((synced + 1))

        elif [ "$ahead_gh" -gt 0 ] && [ "$ahead_or" -gt 0 ]; then
            echo -e "${RED}⚠ CONFLICT${NC}: $branch  (diverged — github +$ahead_gh, origin +$ahead_or — manual merge required)"
            conflicts=$((conflicts + 1))

        else
            skipped=$((skipped + 1))
        fi
    fi
done

echo ""
echo "Syncing tags..."
git push origin --tags --quiet
git push github --tags --quiet

echo ""
echo "Done: $synced branch(es) synced, $skipped already up to date, $conflicts conflict(s)."
